"""
pyrad.io.read_data_sun
======================

Functions for reading data used in sun monitoring

.. autosummary::
    :toctree: generated/

    read_sun_hits_multiple_days
    read_sun_hits
    read_sun_retrieval
    read_solar_flux
"""

import datetime
import csv
from warnings import warn

import numpy as np

from pyart.config import get_fillvalue

from .io_aux import get_save_dir, make_filename


def read_sun_hits_multiple_days(cfg, time_ref, nfiles=1):
    """
    Reads sun hits data from multiple file sources

    Parameters
    ----------
    cfg : dict
        dictionary with configuration data to find out the right file
    time_ref : datetime object
        reference time
    nfiles : int
        number of files to read

    Returns
    -------
    date, ray, nrng, rad_el, rad_az, sun_el, sun_az, ph, ph_std, nph, nvalh,
    pv, pv_std, npv, nvalv, zdr, zdr_std, nzdr, nvalzdr : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    timeinfo = time_ref - datetime.timedelta(days=nfiles - 1)

    date = list()
    ray = list()
    nrng = list()
    rad_el = list()
    rad_az = list()
    sun_el = list()
    sun_az = list()
    ph = np.ma.array(list())
    ph_std = np.ma.array(list())
    nph = list()
    nvalh = list()
    pv = np.ma.array(list())
    pv_std = np.ma.array(list())
    npv = list()
    nvalv = list()
    zdr = np.ma.array(list())
    zdr_std = np.ma.array(list())
    nzdr = list()
    nvalzdr = list()

    for i in range(nfiles):
        savedir = get_save_dir(
            cfg["basepath"],
            cfg["procname"],
            cfg["dsname"],
            cfg["sun_hits_dir"],
            timeinfo=timeinfo,
            create_dir=False,
        )

        fname = make_filename(
            "info",
            cfg["type"],
            "detected",
            ["csv"],
            timeinfo=timeinfo,
            timeformat="%Y%m%d",
        )

        (
            date_aux,
            ray_aux,
            nrng_aux,
            rad_el_aux,
            rad_az_aux,
            sun_el_aux,
            sun_az_aux,
            ph_aux,
            ph_std_aux,
            nph_aux,
            nvalh_aux,
            pv_aux,
            pv_std_aux,
            npv_aux,
            nvalv_aux,
            zdr_aux,
            zdr_std_aux,
            nzdr_aux,
            nvalzdr_aux,
        ) = read_sun_hits(savedir + fname[0])

        if date_aux is None:
            return (
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )

        date += date_aux
        ray += list(ray_aux)
        nrng += list(nrng_aux)
        rad_el += list(rad_el_aux)
        rad_az += list(rad_az_aux)
        sun_el += list(sun_el_aux)
        sun_az += list(sun_az_aux)
        ph = np.ma.append(ph, ph_aux)
        ph_std = np.ma.append(ph_std, ph_std_aux)
        nph += list(nph_aux)
        nvalh += list(nvalh_aux)
        pv = np.ma.append(pv, pv_aux)
        pv_std = np.ma.append(pv_std, pv_std_aux)
        npv += list(npv_aux)
        nvalv += list(nvalv_aux)
        zdr = np.ma.append(zdr, zdr_aux)
        zdr_std = np.ma.append(zdr_std, zdr_std_aux)
        nzdr += list(nzdr_aux)
        nvalzdr += list(nvalzdr_aux)

        timeinfo += datetime.timedelta(days=1)

    ray = np.squeeze(np.array(ray))
    nrng = np.squeeze(np.array(nrng))
    rad_el = np.squeeze(np.array(rad_el))
    rad_az = np.squeeze(np.array(rad_az))
    sun_el = np.squeeze(np.array(sun_el))
    sun_az = np.squeeze(np.array(sun_az))
    nph = np.squeeze(np.array(nph))
    nvalh = np.squeeze(np.array(nvalh))
    npv = np.squeeze(np.array(npv))
    nvalv = np.squeeze(np.array(nvalv))
    nzdr = np.squeeze(np.array(nzdr))
    nvalzdr = np.squeeze(np.array(nvalzdr))

    return (
        date,
        ray,
        nrng,
        rad_el,
        rad_az,
        sun_el,
        sun_az,
        ph,
        ph_std,
        nph,
        nvalh,
        pv,
        pv_std,
        npv,
        nvalv,
        zdr,
        zdr_std,
        nzdr,
        nvalzdr,
    )


def read_sun_hits(fname):
    """
    Reads sun hits data contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date, ray, nrng, rad_el, rad_az, sun_el, sun_az, ph, ph_std, nph, nvalh,
    pv, pv_std, npv, nvalv, zdr, zdr_std, nzdr, nvalzdr : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)

            ray = np.empty(nrows, dtype=int)
            nrng = np.empty(nrows, dtype=int)
            rad_el = np.empty(nrows, dtype=float)
            rad_az = np.empty(nrows, dtype=float)
            sun_el = np.empty(nrows, dtype=float)
            sun_az = np.empty(nrows, dtype=float)
            ph = np.ma.empty(nrows, dtype=float)
            ph_std = np.ma.empty(nrows, dtype=float)
            nph = np.empty(nrows, dtype=int)
            nvalh = np.empty(nrows, dtype=int)
            pv = np.ma.empty(nrows, dtype=float)
            pv_std = np.ma.empty(nrows, dtype=float)
            npv = np.empty(nrows, dtype=int)
            nvalv = np.empty(nrows, dtype=int)
            zdr = np.ma.empty(nrows, dtype=float)
            zdr_std = np.ma.empty(nrows, dtype=float)
            nzdr = np.empty(nrows, dtype=int)
            nvalzdr = np.empty(nrows, dtype=int)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            date = list()
            for i, row in enumerate(reader):
                date.append(
                    datetime.datetime.strptime(row["time"], "%Y-%m-%d %H:%M:%S.%f")
                )
                ray[i] = int(row["ray"])
                nrng[i] = int(row["NPrng"])
                rad_el[i] = float(row["rad_el"])
                rad_az[i] = float(row["rad_az"])
                sun_el[i] = float(row["sun_el"])
                sun_az[i] = float(row["sun_az"])
                ph[i] = float(row["dBm_sun_hit"])
                ph_std[i] = float(row["std(dBm_sun_hit)"])
                nph[i] = int(row["NPh"])
                nvalh[i] = int(row["NPhval"])
                pv[i] = float(row["dBmv_sun_hit"])
                pv_std[i] = float(row["std(dBmv_sun_hit)"])
                npv[i] = int(row["NPv"])
                nvalv[i] = int(row["NPvval"])
                zdr[i] = float(row["ZDR_sun_hit"])
                zdr_std[i] = float(row["std(ZDR_sun_hit)"])
                nzdr[i] = int(row["NPzdr"])
                nvalzdr[i] = int(row["NPzdrval"])

            ph = np.ma.masked_values(ph, get_fillvalue())
            ph_std = np.ma.masked_values(ph_std, get_fillvalue())
            pv = np.ma.masked_values(pv, get_fillvalue())
            pv_std = np.ma.masked_values(pv_std, get_fillvalue())
            zdr = np.ma.masked_values(zdr, get_fillvalue())
            zdr_std = np.ma.masked_values(zdr_std, get_fillvalue())

            csvfile.close()

            return (
                date,
                ray,
                nrng,
                rad_el,
                rad_az,
                sun_el,
                sun_az,
                ph,
                ph_std,
                nph,
                nvalh,
                pv,
                pv_std,
                npv,
                nvalv,
                zdr,
                zdr_std,
                nzdr,
                nvalzdr,
            )

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return (
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )


def read_sun_retrieval(fname):
    """
    Reads sun retrieval data contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    first_hit_time, last_hit_time, nhits_h, el_width_h, az_width_h, el_bias_h,
    az_bias_h, dBm_sun_est, std_dBm_sun_est, sf_h,
    nhits_v, el_width_v, az_width_v, el_bias_v, az_bias_v, dBmv_sun_est,
    std_dBmv_sun_est, sf_v,
    nhits_zdr, zdr_sun_est, std_zdr_sun_est,
    sf_ref, ref_time : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)

            nhits_h = np.empty(nrows, dtype=int)
            el_width_h = np.ma.empty(nrows, dtype=float)
            az_width_h = np.ma.empty(nrows, dtype=float)
            el_bias_h = np.ma.empty(nrows, dtype=float)
            az_bias_h = np.ma.empty(nrows, dtype=float)
            dBm_sun_est = np.ma.empty(nrows, dtype=float)
            std_dBm_sun_est = np.ma.empty(nrows, dtype=float)
            sf_h = np.ma.empty(nrows, dtype=float)

            nhits_v = np.empty(nrows, dtype=int)
            el_width_v = np.ma.empty(nrows, dtype=float)
            az_width_v = np.ma.empty(nrows, dtype=float)
            el_bias_v = np.ma.empty(nrows, dtype=float)
            az_bias_v = np.ma.empty(nrows, dtype=float)
            dBmv_sun_est = np.ma.empty(nrows, dtype=float)
            std_dBmv_sun_est = np.ma.empty(nrows, dtype=float)
            sf_v = np.ma.empty(nrows, dtype=float)

            nhits_zdr = np.empty(nrows, dtype=int)
            zdr_sun_est = np.ma.empty(nrows, dtype=float)
            std_zdr_sun_est = np.ma.empty(nrows, dtype=float)

            sf_ref = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))

            first_hit_time = list()
            last_hit_time = list()
            ref_time = list()
            for i, row in enumerate(reader):
                first_hit_time.append(
                    datetime.datetime.strptime(row["first_hit_time"], "%Y%m%d%H%M%S")
                )
                last_hit_time.append(
                    datetime.datetime.strptime(row["last_hit_time"], "%Y%m%d%H%M%S")
                )

                nhits_h[i] = int(row["nhits_h"])
                el_width_h[i] = float(row["el_width_h"])
                az_width_h[i] = float(row["az_width_h"])
                el_bias_h[i] = float(row["el_bias_h"])
                az_bias_h[i] = float(row["az_bias_h"])
                dBm_sun_est[i] = float(row["dBm_sun_est"])
                std_dBm_sun_est[i] = float(row["std(dBm_sun_est)"])
                sf_h[i] = float(row["sf_h"])

                nhits_v[i] = int(row["nhits_v"])
                el_width_v[i] = float(row["el_width_v"])
                az_width_v[i] = float(row["az_width_v"])
                el_bias_v[i] = float(row["el_bias_v"])
                az_bias_v[i] = float(row["az_bias_v"])
                dBmv_sun_est[i] = float(row["dBmv_sun_est"])
                std_dBmv_sun_est[i] = float(row["std(dBmv_sun_est)"])
                sf_v[i] = float(row["sf_v"])

                nhits_zdr[i] = int(row["nhits_zdr"])
                zdr_sun_est[i] = float(row["ZDR_sun_est"])
                std_zdr_sun_est[i] = float(row["std(ZDR_sun_est)"])

                sf_ref[i] = float(row["sf_ref"])
                if row["ref_time"] == "None":
                    ref_time.append(None)
                else:
                    ref_time.append(
                        datetime.datetime.strptime(row["ref_time"], "%Y%m%d%H%M%S")
                    )

            el_width_h = np.ma.masked_values(el_width_h, get_fillvalue())
            az_width_h = np.ma.masked_values(az_width_h, get_fillvalue())
            el_bias_h = np.ma.masked_values(el_bias_h, get_fillvalue())
            az_bias_h = np.ma.masked_values(az_bias_h, get_fillvalue())
            dBm_sun_est = np.ma.masked_values(dBm_sun_est, get_fillvalue())
            std_dBm_sun_est = np.ma.masked_values(std_dBm_sun_est, get_fillvalue())
            sf_h = np.ma.masked_values(sf_h, get_fillvalue())

            el_width_v = np.ma.masked_values(el_width_v, get_fillvalue())
            az_width_v = np.ma.masked_values(az_width_v, get_fillvalue())
            el_bias_v = np.ma.masked_values(el_bias_v, get_fillvalue())
            az_bias_v = np.ma.masked_values(az_bias_v, get_fillvalue())
            dBmv_sun_est = np.ma.masked_values(dBmv_sun_est, get_fillvalue())
            std_dBmv_sun_est = np.ma.masked_values(std_dBmv_sun_est, get_fillvalue())
            sf_v = np.ma.masked_values(sf_v, get_fillvalue())

            zdr_sun_est = np.ma.masked_values(zdr_sun_est, get_fillvalue())
            std_zdr_sun_est = np.ma.masked_values(std_zdr_sun_est, get_fillvalue())

            sf_ref = np.ma.masked_values(sf_ref, get_fillvalue())

            csvfile.close()

            return (
                first_hit_time,
                last_hit_time,
                nhits_h,
                el_width_h,
                az_width_h,
                el_bias_h,
                az_bias_h,
                dBm_sun_est,
                std_dBm_sun_est,
                sf_h,
                nhits_v,
                el_width_v,
                az_width_v,
                el_bias_v,
                az_bias_v,
                dBmv_sun_est,
                std_dBmv_sun_est,
                sf_v,
                nhits_zdr,
                zdr_sun_est,
                std_zdr_sun_est,
                sf_ref,
                ref_time,
            )

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return (
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )


def read_solar_flux(fname):
    """
    Reads solar flux data from the DRAO observatory in Canada

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    flux_datetime : datetime array
        the date and time of the solar flux retrievals
    flux_value : array
        the observed solar flux

    """
    try:
        with open(fname, "r", newline="") as txtfile:
            # skip the first two lines
            for i in range(2):
                next(txtfile)

            # first count the lines
            reader = csv.DictReader(
                txtfile,
                delimiter=" ",
                skipinitialspace=True,
                fieldnames=[
                    "fluxdate",
                    "fluxtime",
                    "fluxjulian",
                    "fluxcarrington",
                    "fluxobsflux",
                    "fluxadjflux",
                    "fluxursi",
                ],
            )
            nrows = sum(1 for row in reader)

            flux_value = np.empty(nrows, dtype=float)

            # now read the data
            txtfile.seek(0)

            # skip the first two lines
            for i in range(2):
                next(txtfile)

            reader = csv.DictReader(
                txtfile,
                delimiter=" ",
                skipinitialspace=True,
                fieldnames=[
                    "fluxdate",
                    "fluxtime",
                    "fluxjulian",
                    "fluxcarrington",
                    "fluxobsflux",
                    "fluxadjflux",
                    "fluxursi",
                ],
            )

            flux_datetime = list()
            for i, row in enumerate(reader):
                flux_datetime.append(
                    datetime.datetime.strptime(
                        row["fluxdate"] + row["fluxtime"], "%Y%m%d%H%M%S"
                    )
                )
                flux_value[i] = float(row["fluxobsflux"])

            txtfile.close()

            return flux_datetime, flux_value

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None
