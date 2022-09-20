"""
pyrad.io.trajectory
===================
Trajectory class implementation for reading trajectory file.
Converting to different coordinate systems.
.. autosummary::
    :toctree: generated/
    Trajectory
    _Radar_Trajectory
"""

import sys
import re
import datetime
import locale
import pandas as pd
from warnings import warn
from copy import deepcopy

import numpy as np

import pyart

from ..io.read_data_sensor import read_lightning, read_trt_traj_data
from ..io.read_data_sensor import read_trt_thundertracking_traj_data


class Trajectory(object):
    """
    A class for reading and handling trajectory data from a file.
    Attributes
    ----------
    filename : str
        Path and name of the trajectory definition file
    starttime : datetime
        Start time of trajectory processing.
    endtime : datetime
        End time of trajectory processing.
    trajtype : str
    Type of trajectory. Can be 'plane' or 'lightning'
    time_vector : Array of datetime objects
        Array containing the trajectory time samples
    wgs84_lat_deg : Array of floats
        WGS84 latitude samples in radian
    wgs84_lon_deg : Array of floats
        WGS84 longitude samples in radian
    wgs84_alt_m : Array of floats
        WGS84 altitude samples in m
    nsamples : int
    Number of samples in the trajectory
    _swiss_grid_done : Bool
        Indicates that convertion to Swiss coordinates has been performed
    swiss_chy, swiss_chx, swiss_chh : Array of floats
        Swiss coordinates in m
    radar_list : list
          List of radars for which trajectories are going to be computed
    flashnr : int
          For 'lightning' only. Number of flash for which trajectory data
          is going to be computed. If 0 all all flashes are going to be
          considered.
    time_in_flash : array of floats
          For 'lightning' only. Time within flash (sec)
    flashnr_vec : array of ints
          For 'lightning' only. Flash number of each data sample
    dBm : array of floats
          For 'lightning' only. Lightning power (dBm)
    Methods:
    --------
    add_radar : Add a radar
    calculate_velocities : Computes velocities
    get_start_time : Return time of first trajectory sample
    get_end_time : Return time of last trajectory sample
    get_samples_in_period : Get indices of samples within period
    _convert_traj_to_swissgrid : convert data from WGS84 to Swiss coordinates
    _read_traj : Read plane trajectory from file
    _read_traj_flores : Read FLORAKO plane trajectory from file
    _read_traj_trt : Read TRT trajectory from file
    _read_traj_lightning : Read lightning trajectory from file
    _get_total_seconds : Get the total time of the trajectory in seconds
    """

    def __init__(self, filename, starttime=None, endtime=None,
                 trajtype='plane', flashnr=0):
        """
        Initalize the object.
        Parameters
        ----------
        filename : str
            Filename containing the trajectory samples.
        starttime : datetime
            Start time of trajectory processing. If not given, use
            the time of the first trajectory sample.
        endtime : datetime
            End time of trajectory processing. If not given, use
            the time of the last trajectory sample.
        trajtype : str
            type of trajectory. Can be plane or lightning
        flashnr : int
            If type of trajectory is lightning, the flash number to check the
            trajectory. 0 means all flash numbers included
        """

        self.filename = filename
        self.starttime = starttime
        self.endtime = endtime
        self.trajtype = trajtype

        self.time_vector = np.array([], dtype=datetime.datetime)
        self.wgs84_lat_deg = np.array([], dtype=float)
        self.wgs84_lon_deg = np.array([], dtype=float)
        self.wgs84_alt_m = np.array([], dtype=float)
        self.nsamples = None

        self._swiss_grid_done = False
        self.swiss_chy = None
        self.swiss_chx = None
        self.swiss_chh = None

        self.radar_list = []

        if self.trajtype == 'lightning':
            self.flashnr = flashnr
            self.time_in_flash = np.array([], dtype=datetime.datetime)
            self.dBm = np.array([], dtype=float)
            self.flashnr_vec = np.array([], dtype=float)
        elif self.trajtype == 'trt':
            self.cell_contour = np.array([])

        try:
            if self.trajtype == 'lightning':
                self._read_traj_lightning(flashnr)
            elif self.trajtype == 'trt':
                self._read_traj_trt()
            elif self.trajtype == 'plane_flores':
                self._read_traj_flores()
            else:
                self._read_traj()

        except Exception as ee:
            warn(str(ee))
            raise Exception(
                "ERROR: Could not load trajectory data from file '" +
                filename+"' into Trajectory object")

    def add_radar(self, radar):
        """
        Add the coordinates (WGS84 longitude, latitude and non WGS84 altitude)
        of a radar to the radar_list.
        Parameters
        ----------
        radar : pyart radar object
            containing the radar coordinates
        Return
        ------
        Radar object
        """

        # Check if radar location is already in the radar list
        for rad in self.radar_list:
            if (rad.location_is_equal(radar.latitude['data'][0],
                                      radar.longitude['data'][0],
                                      radar.altitude['data'][0])):
                # warn("WARNING: Tried to add the same radar twice to the"
                #      " radar list")
                return rad

        rad = _Radar_Trajectory(radar.latitude['data'][0],
                                radar.longitude['data'][0],
                                radar.altitude['data'][0])
        self.radar_list.append(rad)

        # Convert trajectory WGS84 points to polara radar coordinates
        rad.convert_radpos_to_swissgrid()

        self._convert_traj_to_swissgrid()

        # Note: Earth curvature not considered yet!

        (rvec, azvec, elvec) = pyart.core.cartesian_to_antenna(
            self.swiss_chy-rad.ch_y, self.swiss_chx-rad.ch_x,
            self.swiss_chh-rad.ch_alt)

        rad.assign_trajectory(elvec, azvec, rvec)

        return rad

    def calculate_velocities(self, radar):
        """
        Calculate velocities.
        """

        if not radar.traj_assigned:
            raise Exception("ERROR: No trajectory assigned to radar object")

        dt_secs = np.vectorize(self._get_total_seconds)
        dt = dt_secs(self.time_vector[2:] - self.time_vector[:-2])

        v_r = np.empty(self.nsamples, dtype=float)
        v_r[0] = v_r[-1] = np.nan
        v_r[1:-1] = (radar.range_vec[2:] - radar.range_vec[:-2]) / dt

        v_az = np.empty(self.nsamples, dtype=float)
        v_az[0] = v_az[-1] = np.nan
        daz = radar.azimuth_vec[2:] - radar.azimuth_vec[:-2]
        daz[daz > 180.] -= 360.
        daz[daz < -180.] += 360.
        v_az[1:-1] = daz / dt

        v_el = np.empty(self.nsamples, dtype=float)
        v_el[0] = v_el[-1] = np.nan
        v_el[1:-1] = (radar.elevation_vec[2:] - radar.elevation_vec[:-2]) / dt

        v_abs = np.empty(self.nsamples, dtype=float)
        v_abs[0] = v_abs[-1] = np.nan
        self._convert_traj_to_swissgrid()
        dx = self.swiss_chy[2:] - self.swiss_chy[:-2]
        dy = self.swiss_chx[2:] - self.swiss_chx[:-2]
        dz = self.swiss_chh[2:] - self.swiss_chh[:-2]
        v_abs[1:-1] = np.sqrt(dx**2 + dy**2 + dz**2) / dt

        radar.assign_velocity_vecs(v_abs, v_r, v_el, v_az)

    def get_samples_in_period(self, start=None, end=None):
        """"
        Get indices of samples of the trajectory within given time
        period.
        """

        if ((start is None) and (end is None)):
            raise Exception("ERROR: Either start or end must be defined")
        elif start is None:
            return np.where(self.time_vector < end)
        elif end is None:
            return np.where(self.time_vector >= start)
        else:
            return np.where((self.time_vector >= start) &
                            (self.time_vector < end))

    def get_start_time(self):
        """
        Get time of first trajectory sample.
        Return
        ------
        datetime object
        """

        return self.time_vector[0]

    def get_end_time(self):
        """
        Get time of last trajectory sample.
        Return
        ------
        datetime object
        """

        return self.time_vector[-1]

    def _convert_traj_to_swissgrid(self):
        """
        Convert trajectory samples from WGS84 to Swiss CH1903 coordinates
        """

        if self._swiss_grid_done:
            return

        (self.swiss_chy, self.swiss_chx, self.swiss_chh) = \
            pyart.core.wgs84_to_swissCH1903(
                self.wgs84_lon_deg, self.wgs84_lat_deg, self.wgs84_alt_m)

        self._swiss_grid_done = True

    def _read_traj(self):
        """
        Read trajectory from file
        File format
        -----------
        Comment symbol: '#'
        Columns:
        1. Day (UTC) format: DD-MMM-YYYY
        2. Seconds since midnight (UTC)
        3. WGS84 latitude [radians]
        4. WGS84 longitude [radians]
        5. WGS84 altitude [m]
        """

        # check if the file can be read
        try:
            tfile = open(self.filename, "r")
        except:
            raise Exception("ERROR: Could not find|open trajectory file '" +
                            self.filename+"'")

        repat = re.compile("(\d+\-[A-Za-z]+\-\d+)\s+([\d\.]+)\s+"
                           "([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)")

        try:
            loc_set = False
            loc = locale.getlocale()  # get current locale
            if loc[0] != 'en_US':
                try:
                    locale.setlocale(locale.LC_ALL, ('en_US', 'UTF-8'))
                except Exception as ee:
                    raise Exception("ERROR: Cannot set local 'en_US': %s"
                                    % str(ee))
                loc_set = True

            recording_started = True
            if self.starttime is not None:
                recording_started = False
            recording_check_stop = False
            if self.endtime is not None:
                recording_check_stop = True

            for line in tfile:
                line = line.strip()

                if not line:
                    continue

                # ignore comments
                if line.startswith('#'):
                    continue

                line = line.partition('#')[0]  # Remove comments
                line = line.strip()

                mm = repat.match(line)
                if not mm:
                    print("WARNING: Format error in trajectory file '%s'"
                          " on line '%s'" % (self.filename, line),
                          file=sys.stderr)
                    continue

                # Get time stamp
                try:
                    sday = datetime.datetime.strptime(mm.group(1), "%d-%b-%Y")
                except Exception as ee:
                    print(datetime.datetime.utcnow().strftime("%d-%b-%Y"))
                    raise Exception("ERROR: Format error in traj file '%s' "
                                    "on line '%s' (%s)"
                                    % (self.filename, line, str(ee)))

                sday += datetime.timedelta(seconds=float(mm.group(2)))

                if not recording_started:
                    if sday < self.starttime:
                        continue
                    else:
                        recording_started = True

                if recording_check_stop:
                    if sday > self.endtime:
                        break

                self.time_vector = np.append(self.time_vector, [sday])

                self.wgs84_lat_deg = np.append(
                    self.wgs84_lat_deg, [float(mm.group(3)) * 180. / np.pi])
                self.wgs84_lon_deg = np.append(
                    self.wgs84_lon_deg, [float(mm.group(4)) * 180. / np.pi])
                self.wgs84_alt_m = np.append(
                    self.wgs84_alt_m, [float(mm.group(5))])
        except:
            raise
        finally:
            tfile.close()
            if loc_set:
                locale.setlocale(locale.LC_ALL, loc)  # restore saved locale

        self.nsamples = len(self.time_vector)
        
    def _read_traj_flores(self):
        """
        Read trajectory from FLORAKO file
        File format
        -----------
        Column/Variable Contents, Units and Description:
        01: UTCDate      Year-Month-Day                 Date of epoch or feature (UTC time)
        02: UTCTime      HH:MM:SS.SS                    Time of epoch or feature - Receiver time frame (UTC)
        03: Latitude     Decimal Degrees (signed)       North/South Geographic coordinate
        04: Longitude    Decimal Degrees (signed)       East/West Geographic coordinate
        05: H-MSL        Metres                         Height above the geoid
        06: H-Ell        Metres                         Height above the current ellipsoid
        07: Undulation   Metres                         Height of the ellipsoid above the geoid
        08: PDOP                                        Position Dilution of Precision, which is a measure of X, Y, Z position geometry
        09: HDOP                                        Horizontal Position Dilution of Precision, which is a measure of X, Y position geometry
        10: VDOP                                        Vertical Position Dilution of Precision, which is a measure of Z position geometry
        11: DopRms                                      Root mean square of all L1 Doppler measurement residuals
        """
          # check if the file can be read
        try:
            tfile = open(self.filename, "r")
        except:
            raise Exception("ERROR: Could not find|open trajectory file '" +
                            self.filename+"'")
        tfile.close()
        
        try:
            data = pd.read_csv(self.filename, skiprows=28, 
                               delim_whitespace=True)
            data = data.drop(0)
            times = [datetime.datetime.strptime(d+'_'+h+'0000',
                                                '%Y-%m-%d_%H:%M:%S.%f') for 
                                                d, h in zip(data['UTCDate'],
                                                            data['UTCTime'])]
                       
            self.time_vector = np.array(times)
            self.wgs84_lat_deg = np.array(pd.to_numeric(data['Latitude']))
            self.wgs84_lon_deg = np.array(pd.to_numeric(data['Longitude']))
            self.wgs84_alt_m = np.array(pd.to_numeric(data['H-MSL']))
                
        except:
            raise
            
        self.nsamples = len(self.time_vector)
        
    def _read_traj_lightning(self, flashnr=0):
        """
        Read trajectory from lightning file
        File format
        -----------
        Columns:
        1. flashnr (0=noise)
        2. Seconds since midnight (UTC)
        2. Time in flash (s)
        3. WGS84 latitude [deg]
        4. WGS84 longitude [deg]
        5. WGS84 altitude [m]
        6. Power [dBm]
        Parameters
        ----------
        flashnr : int
            the flash number to keep. If 0 data from all flashes will be kept
        """
        flashnr_vec, time, time_in_flash, lat, lon, alt, dBm = read_lightning(
            self.filename)

        if flashnr_vec is None:
            raise Exception("ERROR: Could not find|open trajectory file '" +
                            self.filename+"'")

        recording_started = True
        if self.starttime is not None:
            recording_started = False
        recording_check_stop = False
        if self.endtime is not None:
            recording_check_stop = True

        if flashnr > 0:
            flashnr_vec_aux = deepcopy(flashnr_vec)
            flashnr_vec = flashnr_vec[flashnr_vec_aux == flashnr]
            time = time[flashnr_vec_aux == flashnr]
            time_in_flash = time_in_flash[flashnr_vec_aux == flashnr]
            lat = lat[flashnr_vec_aux == flashnr]
            lon = lon[flashnr_vec_aux == flashnr]
            alt = alt[flashnr_vec_aux == flashnr]
            dBm = dBm[flashnr_vec_aux == flashnr]

        for i, dBm_val in enumerate(dBm):

            if not recording_started:
                if time[i] < self.starttime:
                    continue
                else:
                    recording_started = True

            if recording_check_stop:
                if time[i] > self.endtime:
                    break

            self.flashnr_vec = np.append(self.flashnr_vec, [flashnr_vec[i]])
            self.time_vector = np.append(self.time_vector, [time[i]])
            self.time_in_flash = np.append(
                self.time_in_flash, [time_in_flash[i]])

            self.wgs84_lat_deg = np.append(self.wgs84_lat_deg, [lat[i]])
            self.wgs84_lon_deg = np.append(self.wgs84_lon_deg, [lon[i]])
            self.wgs84_alt_m = np.append(self.wgs84_alt_m, [alt[i]])

            self.dBm = np.append(self.dBm, [dBm_val])

        self.nsamples = len(self.time_vector)

    def _read_traj_trt(self):
        """
        Read trajectory from TRT file
        File format
        -----------
        Columns:
        1. traj_ID
        2. yyyymmddHHMM (UTC)
        2. lon (deg)
        3. lat [deg]
        4. ell_L
        5. ell_S
        6. ell_or
        7. area
        8. vel_x
        9. vel_y
        10. det
        11. RANKr
        12. CG_n
        13. CG_p
        14. CG
        15. CG_percent_p
        16. ET45
        17. ET45m
        18. ET15
        19. ET15m
        20. VIL
        21. maxH
        22. maxHm
        23. POH
        24. RANK
        25. Dvel_x
        26. Dvel_y
        27. cell_contours
        Parameters
        ----------
        """
        if '_tt.trt' in self.filename:
            (traj_ID, _, yyyymmddHHMM, _, _, _, lon, lat, _, _, _, _, _, _, _,
             _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
             cell_contours) = read_trt_thundertracking_traj_data(
                 self.filename)

            if traj_ID is None:
                raise Exception(
                    "ERROR: Could not find|open trajectory file '" +
                    self.filename+"'")

            valid = np.logical_not(np.ma.getmaskarray(yyyymmddHHMM))
            yyyymmddHHMM = yyyymmddHHMM[valid]
            traj_ID = traj_ID[valid]
            lon = lon[valid]
            lat = lat[valid]
            cell_contours = cell_contours[valid]

            if traj_ID.size == 0:
                raise Exception(
                    "ERROR: No valid data in trajectory file '" +
                    self.filename+"'")

        else:
            (traj_ID, yyyymmddHHMM, lon, lat, _, _, _, _, _, _, _, _, _, _, _,
             _, _, _, _, _, _, _, _, _, _, _, _, cell_contours) = (
                 read_trt_traj_data(self.filename))

            if traj_ID is None:
                raise Exception(
                    "ERROR: Could not find|open trajectory file '" +
                    self.filename+"'")

        recording_started = True
        if self.starttime is not None:
            recording_started = False
        recording_check_stop = False
        if self.endtime is not None:
            recording_check_stop = True

        for i, cell_contour in enumerate(cell_contours):
            if not recording_started:
                if yyyymmddHHMM[i] < self.starttime:
                    continue
                else:
                    recording_started = True

            if recording_check_stop:
                if yyyymmddHHMM[i] > self.endtime:
                    break

            self.time_vector = np.append(self.time_vector, [yyyymmddHHMM[i]])

            self.wgs84_lat_deg = np.append(self.wgs84_lat_deg, [lat[i]])
            self.wgs84_lon_deg = np.append(self.wgs84_lon_deg, [lon[i]])
            self.wgs84_alt_m = np.append(self.wgs84_alt_m, 0.)

            self.cell_contour = np.append(self.cell_contour, [cell_contour])

        self.nsamples = len(self.time_vector)

    def _get_total_seconds(self, x):
        """ Return total seconds of timedelta object"""
        return x.total_seconds()


class _Radar_Trajectory:
    """
    A class for holding the trajectory data assigned to a radar.
    Attributes
    ----------
    latitude : float
       WGS84 radar latitude [deg]
    longitude : float
       WGS84 radar longitude [deg]
    altitude : float
       radar altitude [m] (non WGS84)
    ch_y, ch_x, ch_alt : float
       radar coordinates in swiss CH1903 coordinates
    elevation_vec : float list
       Elevation values of the trajectory samples
    azimuth_vec : float list
       Azimuth values of the trajectory samples
    range_vec : float list
       Range values of the trajectory samples
    v_abs, v_r, v_el, v_az : array-like
       Velocity vectors of the absolute [m/s], radial [m/s], elevation [deg/s]
       and azimuth [deg/s] velocities
    Methods:
    --------
    location_is_equal
    assign_trajectory
    convert_radpos_to_swissgrid
    assign_velocity_vecs
    """

    def __init__(self, lat, lon, alt):
        """
        Initalize the object.
        Parameters
        ----------
        lat, lon , alt : radar location coordinates
        nsamps : number of samples
        """

        # radar position in lat, lon
        self.latitude = lat
        self.longitude = lon
        self.altitude = alt

        # radar position in swiss CH1903 coordinates
        self._swiss_coords_done = False
        self.ch_y = None
        self.ch_x = None
        self.ch_alt = None

        # trajectory in polar radar coordinates
        self.traj_assigned = False
        self.elevation_vec = None
        self.azimuth_vec = None
        self.range_vec = None

        # Velocity vectors
        self._velocity_vecs_assigned = False
        self.v_abs = None
        self.v_r = None
        self.v_el = None
        self.v_az = None

    def location_is_equal(self, lat, lon, alt):
        """
        Check if the given coordinates are the same.
        Parameters
        ----------
        lat, lon , alt : radar location coordinates
        Return
        ------
        True if the radar location is equal, False otherwise
        """

        lat_tol = 0.002  # [deg]
        lon_tol = 0.002  # [deg]
        alt_tol = 2.0    # [m]

        if ((np.abs(self.latitude-lat) > lat_tol) or
                (np.abs(self.longitude-lon) > lon_tol) or
                (np.abs(self.altitude-alt) > alt_tol)):
            return False

        return True

    def assign_trajectory(self, el, az, rr):
        """
        Assign a trajectory to the radar in polar radar
        coordinates.
        Parameters
        ----------
        el, az, rr : array-like
           elevation, azimuth and range vector
        """

        if self.traj_assigned:
            warn("WARNING: Trajectory already assigned")
            return

        self.elevation_vec = el
        self.azimuth_vec = az
        self.range_vec = rr
        self.traj_assigned = True

    def convert_radpos_to_swissgrid(self):
        """
        Convert the radar location (in WGS84 coordinates) to
        swiss CH1903 coordinates.
        """

        if self._swiss_coords_done:
            return

        (self.ch_y, self.ch_x, self.ch_alt) = \
            pyart.core.wgs84_to_swissCH1903(self.longitude,
                                            self.latitude,
                                            self.altitude,
                                            no_altitude_transform=True)
        self._swiss_coords_done = True

    def assign_velocity_vecs(self, v_abs, v_r, v_el, v_az):
        """
        Assign velocity vectors to the radar.
        """

        if self._velocity_vecs_assigned:
            warn("WARNING: Trajectory velocity vectors already assigned")
            return

        self.v_abs = v_abs
        self.v_r = v_r
        self.v_az = v_az
        self.v_el = v_el
        self._velocity_vecs_assigned = True
