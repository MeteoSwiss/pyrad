/*
  *************************************************************

  Read rainbow files

  *************************************************************

  Filename:         rainbow_read_raw.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-19
  Last update:      2015-01-26 (fvj: added reading of sensor info)

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Read data and parameters from rainbow files.

  *************************************************************
  =============================================================
 */
/*
* This is a simple program to read raw Selex 50DX files (Rainbow format)
* It accepts all 3 types of scans (volume, elevation or azimuth), for all
* variables available on the system.
*
* It simply loops through the header and parses the values to some structs,
* then reads the BLOB (binary large object as they call) and write to a file.
*
* It reads only Qt-compressed Rainbow files (which is the standard)
*
* The output file is structured as follows:
*
* For AZI/ELE files:
*
* Number of slices     NSLIC (unsigned short integer, 2 bytes)
*  - Always 1 as there is only one slice for each azi/ele file
* Number of angles     NANGL (unsigned short integer, 2 bytes)
* Azimuths matrix      ANGL[NANGL] (float, 4 bytes)
* Number of bins       NBINS (unsigned short integer, 2 bytes)
* Data matrix          DATA[NANGL*NBINS] (float, 4 bytes)
*
* For VOL files:
*
* Number of slices     NSLIC (unsigned short integer, 2 bytes)
*  - Number of slices = number of elevations
*
* Number of angles (slice 1)    NANGL (unsigned short integer, 2 bytes)
* Azimuths matrix (slice 1)     ANGL[NANGL] (float, 4 bytes)
* Number of bins (slice 1)      NBINS (unsigned short integer, 2 bytes)
* Data matrix (slice 1)         DATA[NANGL*NBINS] (float, 4 bytes)
* Number of angles (slice 2)    NANGL (unsigned short integer, 2 bytes)
* Azimuths matrix (slice 2)     ANGL[NANGL] (float, 4 bytes)
* Number of bins (slice 2)      NBINS (unsigned short integer, 2 bytes)
* Data matrix (slice 2)         DATA[NANGL*NBINS] (float, 4 bytes)
* ...
* Number of angles (slice N)    NANGL (unsigned short integer, 2 bytes)
* Azimuths matrix (slice N)     ANGL[NANGL] (float, 4 bytes)
* Number of bins (slice N)      NBINS (unsigned short integer, 2 bytes)
* Data matrix (slice N)         DATA[NANGL*NBINS] (float, 4 bytes)
*
*
* Thiago Biscaro - thiago.biscaro@cptec.inpe.br
*
* */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "endianness.h"
#include "vars_def.h"
#include "qUncompress.h"

#ifdef __cplusplus
extern "C" {
#endif

/* forward definition */
static void return_value(char *source, char *result);
static void fill_dyn_info(char *line, struct params_slice *slice);
static void return_info(char *line, char *pattern, char *result);
static int readFileHeader(char *filename);
static int getBlobData(char *fname, int myblobid, float *data, long nvals);

/* static variables */
static struct blob_ray_info blob_info[MAX_BLOBS];
static struct params_slice p_slice[MAX_SLICES];
static struct params_radar p_radar;
static struct params_sensor p_sensor;
static char scan_type = 0;
static fpos_t fileposition;

/* FUNCTIONS TO READ SENSOR INFORMATION */

/**
 * Get sensor longitude from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return sensor longitude or -99999. if an error occured.
 */
float rainbow_getLongitude(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }
  return p_sensor.lon;
}

/**
 * Get sensor latitude from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return number of slices or -99999. if an error occured.
 */
float rainbow_getLatitude(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }
  return p_sensor.lat;
}

/**
 * Get sensor altitude from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return number of slices or -99999. if an error occured.
 */
float rainbow_getAltitude(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }
  return p_sensor.alt;
}

/**
 * Get sensor wavelength from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return number of slices or -99999. if an error occured.
 */
float rainbow_getWavelength(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }
  return p_sensor.wavelen;
}

/**
 * Get sensor beamwidth from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return number of slices or -99999. if an error occured.
 */
float rainbow_getBeamwidth(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }
  return p_sensor.beamwidth;
}

/**
 * Get sensor type
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  str       pointer to sensor type string
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSensorType(char *fname, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  //printf("sensor type: %s\r\n", p_sensor.type);
  *str = p_sensor.type;
  return 0;
}

/**
 * Get sensor id
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  str       pointer to sensor ID string
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSensorID(char *fname, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  *str = p_sensor.id;
  return 0;
}

/**
 * Get sensor name
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  str       pointer to sensor name string
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSensorName(char *fname, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  *str = p_sensor.name;
  return 0;
}

/* END OF FUNCTIONS TO READ SENSOR INFORMATION */


/**
 * Get horizontal channel nosie power [dBZ at 1 km]
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 *
 * @return noise power on success, -99999. otherwise
 */
float rainbow_getNoisePowerZh(char *fname, int slicenum)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -99999.;
  }

  return p_slice[slicenum].noiseh;
}

/**
 * Get vertical channel nosie power [dBZ at 1 km]
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 *
 * @return noise power on success, -99999. otherwise
 */
float rainbow_getNoisePowerZv(char *fname, int slicenum)
{
  if (readFileHeader(fname) < 0) {
    return -99999.;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -99999.;
  }

  return p_slice[slicenum].noisev;
}

/**
 * Get horizontal radar calibration constants
 *
 * @param fname     Name of file to read (.vol, .azi or .ele)
 * @param slicenum	 slice number
 * @param radarconstants  Radar constants array
 * @param len Length of radar constants array
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getRadarConstanth(char *fname, int slicenum, float *radarconstants, int len)
{
  int k;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  for (k=0; k<len; k++) {
    radarconstants[k] = p_slice[slicenum].radconsth[k];
  }
  return 0;
}

/**
 * Get vertical radar calibration constants
 *
 * @param fname     Name of file to read (.vol, .azi or .ele)
 * @param slicenum	 slice number
 * @param radarconstants  Radar constants array
 * @param len Length of radar constants array
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getRadarConstantv(char *fname, int slicenum, float *radarconstants, int len)
{
  int k;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  for (k=0; k<len; k++) {
    radarconstants[k] = p_slice[slicenum].radconstv[k];
  }
  return 0;
}


/**
 * Get pulse width index from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return pulse width index or -1 if an error occured.
 */
int rainbow_getPulseWidthIndex(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  return p_radar.pw_index;
}



/**
 * Get number of slices from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return number of slices or -1 if an error occured.
 */
int rainbow_getNumberOfSlices(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  return p_radar.numele;
}

/**
 * Get scan name
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  str       pointer to timestr
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getScanName(char *fname, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  *str = p_radar.scanname;
  return 0;
}

/**
 * Get range resolution in km
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return Range resolution [km], -1 on error
 */
float rainbow_getRangeResolution(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  return p_radar.r_s.rangestep;
}

/**
 * Get range sampling
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return Range Sampling, -1 on error
 */
int rainbow_getRangeSampling(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  return p_radar.r_s.rangesamp;
}

/**
 * Get angle step
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 *
 * @return Angle Step [Deg], -1 on error
 */
float rainbow_getAngleStep(char *fname)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }
  return p_radar.r_s.anglestep;
}

/**
 * Get slice angles (= elevation angles for a volume scan) from rainbow file
 *
 * @param  fname  Name of file to read (.vol, .azi or .ele)
 * @param  sliceangles Elevation angles
 * @param  len    Length of sliceangles array
 *
 * @return 0 on sucess, -1 otherwise
 */
int rainbow_getSliceAngles(char *fname, float *sliceangles, int len)
{
  int k;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  if (len != p_radar.numele) {
    printf("%s:%s: Warning: Length of slicearray does not correspond to number of slices.\r\n",
	   __FILE__, __func__);
  }
  for (k=0; k<len; k++) {
    sliceangles[k] = p_slice[k].angle;
  }

  return 0;
}

/**
 * Get slice date
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  str       pointer to datestr
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSliceDate(char *fname, int slicenum, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  *str = p_slice[slicenum].date;
  return 0;
}

/**
 * Get slice time
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  str       pointer to timestr
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSliceTime(char *fname, int slicenum, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  *str = p_slice[slicenum].time;
  return 0;
}

/**
 * Get slice datatype
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  str       pointer to typestr
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSliceDatatype(char *fname, int slicenum, char **str)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  *str = p_slice[slicenum].typestr;
  return 0;
}

/**
 * Get fixed angle of slice
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  pointer to fixed angle
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getFixedAngle(char *fname, int slicenum, float *fixedangle)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  *fixedangle = p_slice[slicenum].angle;
  return 0;
}

/**
 * Get antenna speed during slice recording.
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  pointer to antspeed
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getAntSpeed(char *fname, int slicenum, float *antspeed)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n",
           __FILE__, __func__, slicenum);
    return -1;
  }

  *antspeed = p_slice[slicenum].r_s.antspeed;
  return 0;
}

/**
 * Get number of angles of a slice of a rainbow file
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 *
 * @return number of angles or -1 if an error occured.
 */
int rainbow_getNumberOfAngles(char *fname, int slicenum)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  return p_slice[slicenum].rays;
}

/**
 * Get number of range bins of a slice of a rainbow file
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 *
 * @return number of bins or -1 if an error occured.
 */
int rainbow_getNumberOfRangeBins(char *fname, int slicenum)
{
  if (readFileHeader(fname) < 0) {
    return -1;
  }

  if (slicenum >= p_radar.numele) {
    printf("%s:%s: ERROR: Slice \'%d\' does not exist\r\n", __FILE__, __func__, slicenum);
    return -1;
  }

  return p_slice[slicenum].bins;;
}

/**
 * Get start angles of a slice
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  data      Pointer to data to store
 * @param  nangles   Number of angles
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getStartAngles(char *fname, int slicenum, float *data, long nangles)
{
  int blobid;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  blobid = p_slice[slicenum].blobid_startangle;
  if (blobid < 0) {
    printf("%s %s ERROR: No blobid for startangle of slice %d found.\r\n", __FILE__, __func__, slicenum);
    return -1;
  }
  return getBlobData(fname, blobid, data, nangles);
}

/**
 * Get stop angles of a slice
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  data      Pointer to data to store
 * @param  nangles   Number of angles
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getStopAngles(char *fname, int slicenum, float *data, long nangles)
{
  int blobid;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  blobid = p_slice[slicenum].blobid_stopangle;
  if (blobid < 0) {
    printf("%s %s ERROR: No blobid for stopangle of slice %d found.\r\n", __FILE__, __func__, slicenum);
    return -1;
  }
  return getBlobData(fname, blobid, data, nangles);
}

/**
 * Get data of a slice
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  data      Pointer to data to store
 * @param  nvals     Number of data values
 *
 * @return 0 on success, -1 otherwise
 */
int rainbow_getSliceData(char *fname, int slicenum, float *data, long nvals)
{
  int blobid;

  if (readFileHeader(fname) < 0) {
    return -1;
  }
  blobid = p_slice[slicenum].blobid;
  if (blobid < 0) {
    printf("%s %s  ERROR: No blobid for data of slice %d found.\r\n", __FILE__, __func__, slicenum);
    return -1;
  }
  return getBlobData(fname, blobid, data, nvals);
}

/**
 * Get data of a blob
 *
 * @param  fname     Name of file to read (.vol, .azi or .ele)
 * @param  slicenum  Number of the slice (0=first slice)
 * @param  data      Pointer to data to store
 * @param  nvals     Number of data values
 *
 * @return 0 on success, -1 otherwise
 */
  static int getBlobData(char *fname, int myblobid, float *data, long nvals)
  {
    FILE *fp;
    char line[MAX_LINE_SIZE];
    char line2[MAX_LINE_SIZE];
    char *buf = NULL;
    unsigned short int *buffer_16 = NULL;
    unsigned char *buffer_8 = NULL;
    size_t nbytes = MAX_LINE_SIZE - 1;
    char tmp_str[MAX_LINE_SIZE];
    int i;
    unsigned short int blobid = 0;
    unsigned int size_blob = 0;
    float step = 0;
    float max = 0, min = 0;
    float value_obs = 0;
    size_t len;

    if (NULL == (fp = fopen(fname, "r"))) {
      printf("%s:%s:  Error opening input file %s\r\n", __FILE__, __func__, fname);
      return -1;
    }
    /* Set file position after header was read. */
    fsetpos(fp, &fileposition);

    /* -------------------------------------------------------------------------*/
    /* --------- Data reading                                     --------------*/
    /* -------------------------------------------------------------------------*/
    while (feof(fp) == 0)
      {
        /*reads binary information (and headers)*/
        memset(line, 0, MAX_LINE_SIZE);
        if (fgets(line2, nbytes, fp) == NULL) {
          if (feof(fp) == 1) {
            //printf("%s:%s:  End of file reached\r\n", __FILE__, __func__);
            fclose(fp);
            return 0;
          }
          if (ferror(fp) == 1) {
            printf("%s:%s:  ERROR: error reading file\r\n", __FILE__, __func__);
            fclose(fp);
            return -1;
          }
        }
        strncpy(line, line2, MAX_LINE_SIZE);
        /*Now reads BLOB ID info*/

        if (NULL != strstr(line, "<BLOB "))
          {
            strcpy(tmp_str, "blobid");
            return_info(line, tmp_str, tmp_str);
            blobid = atoi(tmp_str);

            //slice = blob_info[blobid].slice;
            //depth = blob_info[blobid].depth;
            strcpy(tmp_str, "size");
            return_info(line, tmp_str, tmp_str);
            size_blob = atoi(tmp_str);
            strcpy(tmp_str, "compression");
            return_info(line, tmp_str, tmp_str);

            /* Check if is the right slice */
            if (blobid != myblobid) {
              /* */
              fseek(fp, size_blob, SEEK_CUR);
              continue;
            }

            /*reads binary data*/
            buf = (char *) malloc(size_blob);
            if (NULL == buf) {
              fclose(fp);
              printf("%s:%s ERROR: Blob size not found\r\n", __FILE__, __func__);
              return -1;
            }
            len = fread(buf, size_blob, 1, fp);

            //printf("Read data from blobid %d\r\n", blobid);
            if (len)
              {
                if (qUncompress((const unsigned char *)buf, size_blob, &buffer_8, &size_blob) < 0) {
                  printf("%s:%s:  ERROR uncompressing data\r\n", __FILE__, __func__);
                  fclose(fp);
                  return -1;
                }

                if (blob_info[blobid].depth == 16) {
                  /*16 bit data (put on Little Endian order)*/
                  buffer_16 = (unsigned short int *) malloc(size_blob);

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
                  unsigned k;
                  for (k = 0; k < size_blob - 1; k = k + 2)
                    ((unsigned char *)buffer_16)[k] = buffer_8[k+1];
                  for (k = 1; k < size_blob; k = k + 2)
                    ((unsigned char *)buffer_16)[k] = buffer_8[k-1];
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
		  memcpy(buffer_16, buffer_8, size_blob);
#else
#error "Endianness undefined!"
#endif
                }

                if (blob_info[blobid].type == BLOB_STARTANGLE ||
                    blob_info[blobid].type == BLOB_STOPANGLE) {
                  //printf("Reading slice %d (angles)\r\n", blob_info[blobid].slice + 1);

                  /* Check data length */
                  if (nvals != p_slice[blob_info[blobid].slice].rays) {
                    printf("%s:%s:  ERROR: Unexpected number of angles (expected:%ld, found:%d)\r\n",
                           __FILE__, __func__,
                           nvals, p_slice[blob_info[blobid].slice].rays);
                    fclose(fp);
                    return -1;
                  }

                  /*always 16 bits*/
                  step = 360.0/65536.0;
                  for (i = 0; i < nvals; i++) {
                    value_obs = (float) (buffer_16[i] * step);
                    /*
                      from Rainbow File Format Manual:
                      if angle is greater than 225 degrees, subtract 360.
                      For the 50DX the valid range for elevation scans is
                      -2 to 182 degrees. So 358 become -2 degrees.
                    */
                    if ((scan_type == SCAN_RHI) && value_obs > 225)
                      value_obs = value_obs - 360;

                    *(data++) = value_obs;
                  }
                }

                if (blob_info[blobid].type >= BLOB_RAWDATAUZ &&
                    blob_info[blobid].type <= BLOB_RAWDATALDR) {
                  switch (blob_info[blobid].type) {
                  case BLOB_RAWDATAUZ:
                  case BLOB_RAWDATACZ:
                    {
                      min = p_slice[blob_info[blobid].slice].dynz[0];
                      max = p_slice[blob_info[blobid].slice].dynz[1];
                      break;
                    }
                  case BLOB_RAWDATAV:
                    {
                      min = p_slice[blob_info[blobid].slice].dynv[0];
                      max = p_slice[blob_info[blobid].slice].dynv[1];
                      break;
                    }
                  case BLOB_RAWDATAW:
                    {
                      min = p_slice[blob_info[blobid].slice].dynw[0];
                      max = p_slice[blob_info[blobid].slice].dynw[1];
                      break;
                    }
                  case BLOB_RAWDATAKDP:
                  case BLOB_RAWDATAUKDP:
                    {
                      min = p_slice[blob_info[blobid].slice].dynkdp[0];
                      max = p_slice[blob_info[blobid].slice].dynkdp[1];
                      break;
                    }
                  case BLOB_RAWDATAPDP:
                  case BLOB_RAWDATAUPDP:
                    {
                      min = 0;
                      max = 360;
                      break;
                    }
                  case BLOB_RAWDATAZDR:
                    {
                      min = p_slice[blob_info[blobid].slice].dynzdr[0];
                      max = p_slice[blob_info[blobid].slice].dynzdr[1];
                      break;
                    }
                  case BLOB_RAWDATALDR:
                    {
                      min = p_slice[blob_info[blobid].slice].dynldr[0];
                      max = p_slice[blob_info[blobid].slice].dynldr[1];
                      break;
                    }
                  case BLOB_RAWDATARHV:
                  case BLOB_RAWDATASQI:
                    {
                      min = 0;
                      max = 1;
                      break;
                    }
                  }

                  //printf("Reading slice %d (data)\r\n", blob_info[blobid].slice + 1);

                  /* Check data length */
                  if (nvals != p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins) {
                    printf( "%s:%s: ERROR: Unexpected number of data values (expected:%ld, found:%d)\r\n",
                            __FILE__, __func__,
                            nvals, p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins);
                    fclose(fp);
                    return -1;
                  }

                  /*according to Gematronik documentation*/
                  if (blob_info[blobid].depth == 8)
                    step = (max - min)/254;
                  else
                    step = (max - min)/65534;

                  for (i = 0; i < nvals; i++) {
                    value_obs = -9999;
                    switch (blob_info[blobid].depth) {
                    case 16:
                      if (buffer_16[i]) {
                        value_obs = (float) ((buffer_16[i]-1) * step + min);
                      }
                      break;
                    case 8:
                      if (buffer_8[i]) {
                        value_obs = (float) ((buffer_8[i]-1) * step + min);
                      }
                      break;
                    }
                    *(data++) = value_obs;
                  }
                }

                if (blob_info[blobid].type == BLOB_COSMOISO0 ||
                    blob_info[blobid].type == BLOB_COSMOTEMP) {
                  switch (blob_info[blobid].type) {
                  case BLOB_COSMOISO0:
                    {
                      min = 0;
                      max = 254;
                      break;
                    }
                  case BLOB_COSMOTEMP:
                    {
                      min = -87.;
                      max = 40.;
                      break;
                    }
                  }

                  /* Check data length */
                  if (nvals != p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins) {
                    printf( "%s:%s: ERROR: Unexpected number of data values (expected:%ld, found:%d)\r\n",
                            __FILE__, __func__,
                            nvals, p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins);
                    fclose(fp);
                    return -1;
                  }

                  /*according to Gematronik documentation*/
                  if (blob_info[blobid].depth == 8)
                    step = (max - min)/254;
                  else
                    step = (max - min)/65534;

                  for (i = 0; i < nvals; i++) {
                    value_obs = -9999;
                    switch (blob_info[blobid].depth) {
                    case 16:
                      if (buffer_16[i]) {
                        value_obs = (float) ((buffer_16[i]-1) * step + min);
                      }
                      break;
                    case 8:
                      if (buffer_8[i]) {
                        value_obs = (float) ((buffer_8[i]-1) * step + min);
                      }
                      break;
                    }
                    *(data++) = value_obs;
                  }
                }

                if (blob_info[blobid].type == BLOB_COORD) {
                  /* Check data length */
                  if (nvals != p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins) {
                    printf( "%s:%s: ERROR: Unexpected number of data values (expected:%ld, found:%d)\r\n",
                            __FILE__, __func__,
                            nvals, p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins);
                    fclose(fp);
                    return -1;
                  }
                  min = -32767;
                  max = 32767;
                  step=(max-min)/65534;

                  for (i = 0; i < nvals*3; i++) {
                    value_obs = -9999;
                    switch (blob_info[blobid].depth) {
                    case 16:
                      if (buffer_16[i]) {
                        value_obs = (float) (buffer_16[i]-1)*step+min;
                      }
                      break;
                    case 8:
                      if (buffer_8[i]) {
                        value_obs = (float) (buffer_8[i]-1)*step+min;
                      }
                      break;
                    }
                    *(data++) = value_obs;
                  }
                }

                if (blob_info[blobid].type == BLOB_VISIBILITY) {
                  min = 0;
                  max = 127;

                  /* Check data length */
                  if (nvals != p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins) {
                    printf( "%s:%s: ERROR: Unexpected number of data values (expected:%ld, found:%d)\r\n",
                            __FILE__, __func__,
                            nvals, p_slice[blob_info[blobid].slice].rays*p_slice[blob_info[blobid].slice].bins);
                    fclose(fp);
                    return -1;
                  }

                  /*according to Gematronik documentation*/
                  if (blob_info[blobid].depth == 8)
                    step = (max - min)/254;
                  else
                    step = (max - min)/65534;

                  for (i = 0; i < nvals; i++) {
                    value_obs = -9999;
                    switch (blob_info[blobid].depth) {
                    case 16:
                      if (buffer_16[i]) {
                        value_obs = (float) ((buffer_16[i]-1) * step + min);
                      }
                      break;
                    case 8:
                      if (buffer_8[i]) {
                        value_obs = (float) ((buffer_8[i]-1) * step + min);
                      }
                      break;
                    }
                    *(data++) = value_obs;
                  }
                }

                free(buffer_8);

                if (blob_info[blobid].depth == 16)
                  free(buffer_16);
              }
            else {
              printf("%s:%s: ERROR: Unable to read blob data\r\n", __FILE__, __func__);
              fclose(fp);
              return -1;
            }
            free(buf);
          }
      }

    fclose(fp);
    return 0;
  }


static void return_value(char *source, char *result)
{
  int i=0, j=0;

  memset(result, 0, MAX_LINE_SIZE);

  for (i = 0; i < MAX_LINE_SIZE; i++)
    {
      if (source[i] == '>')
	{
	  while (source[++i] != '<')
            {
	      result[j++] = source[i];
            }
	  return;
	}
    }
}

/**
 * Read header of rainbow raw file
 *
 * @param fname  Filename (with path) of rainbow file
 *
 * @return 0 on success, -1 otherwise
 */
static int readFileHeader(char *fname) {

	#define MAXFILENAME 1024
	static char filename[MAXFILENAME] = "";
	static unsigned char fileheader_read = 0;

	FILE *fp;
	char line[MAX_LINE_SIZE];
	char line2[MAX_LINE_SIZE];
	unsigned short int stagger = 0;
	char tmp_str[MAX_LINE_SIZE];
	char * tmp_str2;
	int slice = -1;
	//char begin_pargroup = 0;
	char end_pargroup = 0;
	//char begin_slice = 0;
	//char end_slice = 0;
	int i;
	unsigned short int blobid = 0;
	//struct params_sensor p_sensor;

	/* Read file header if not yet read */
	if ((fileheader_read > 0) && (strcmp(fname, filename) == 0) ) {
		/* header already read */
		return 0;
	}
	fileheader_read = 0;
	filename[0] = '\0';

	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("Error opening input file %s\r\n", fname);
		return -1;
	}

	memset(&p_radar, 0, sizeof(struct params_radar));
	memset(p_slice, 0, MAX_SLICES * sizeof(struct params_slice));
	memset(&p_sensor, 0, sizeof(struct params_sensor));
	memset(blob_info, 0, MAX_BLOBS * sizeof(struct blob_ray_info));
	for (i=0;i<MAX_SLICES;i++) {
		p_slice[i].blobid = -1;
		p_slice[i].blobid_startangle = -1;
		p_slice[i].blobid_stopangle = -1;
	}

	memset(line, 0, MAX_LINE_SIZE);
	while (NULL == strstr(line, "<!-- END XML -->")) {
		memset(line, 0, MAX_LINE_SIZE);

		if (fgets(line2, MAX_LINE_SIZE , fp) == NULL) {
			fclose(fp);
			return -1;
		}
		strncpy(line, line2, MAX_LINE_SIZE);
		if (NULL != strstr(line, "type=\"vol\"")) {
			//printf("Volume found\r\n");
			scan_type = SCAN_VOL;
		}
		if (NULL != strstr(line, "type=\"ele\"")) {
			//printf("Elevation scan found\r\n");
			scan_type = SCAN_RHI;
		}
		if (NULL != strstr(line, "type=\"azi\"")) {
			//printf("Azimuth scan found\r\n");
			scan_type = SCAN_PPI;
		}

		if (NULL != strstr(line, "<scan ")) {
			strcpy(tmp_str, "name");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_radar.scanname, tmp_str, MAX_SCANNAME_LEN);
		}

		if (NULL != strstr(line, "<pargroup ")) {
			//begin_pargroup = 1;
			/*beginning of pargroup*/
		}

		if (NULL != strstr(line, "</pargroup>")) {
			/*end of pargroup*/
			end_pargroup = 1;
			/*copying common info to all slices*/
			for (i = 0; i < MAX_SLICES - 1; i++) {
				memcpy(&(p_slice[i].r_s), &(p_radar.r_s),
					sizeof(struct common_radar_slice));
			}
		}

		if (NULL != strstr(line, "<slice ")) {
			/*beginning of slice*/
			//begin_slice = 1;
			slice++;
			strcpy(tmp_str, "refid");
			return_info(line, tmp_str, tmp_str);
			if ((atoi(tmp_str) != slice))
			{
				printf("Invalid slice index (was looking for %d, found %d)\r\n",
					slice, atoi(tmp_str));
				fclose(fp);
				return -1;
			}
		}
		if (NULL != strstr(line, "</slice>")) {
			/*beginning of slice*/
			//end_slice = 1;
		}
		if (NULL != strstr(line, "<stoprange>")) {
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.stoprange = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.stoprange = (float) atof(tmp_str);
		}

		if (NULL != strstr(line, "<numele>")) {
			return_value(line, tmp_str);
			p_radar.numele = atoi(tmp_str);
			//printf("Found %d slice(s)\r\n", p_radar.numele);
		}
		if (NULL != strstr(line, "<rangestep>")) {
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.rangestep = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.rangestep = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<startrange>")) {
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.startrange = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.startrange = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<highprf>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.highprf = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.highprf = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<lowprf>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.lowprf = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.lowprf = (float) atof(tmp_str);
		}

		if (NULL != strstr(line, "<antspeed>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.antspeed = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.antspeed = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<anglestep>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.anglestep = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.anglestep = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<timesamp>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.timesamp = (float) atof(tmp_str);
			else
				p_slice[slice].r_s.timesamp = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<rangesamp>"))
		{
			return_value(line, tmp_str);
			if (end_pargroup == 0)
				p_radar.r_s.rangesamp = atoi(tmp_str);
			else
				p_slice[slice].r_s.rangesamp = atoi(tmp_str);
		}

		if (NULL != strstr(line, "<posangle>"))
		{
			return_value(line, tmp_str);
			p_slice[slice].angle = (float) atof(tmp_str);
		}

		if (NULL != strstr(line, "<stagger>"))
		{
			if (NULL != strstr(line, "None"))
				stagger = 0;
			if (NULL != strstr(line, "5/4"))
				stagger = 5;
			if (NULL != strstr(line, "4/3"))
				stagger = 4;
			if (NULL != strstr(line, "3/2"))
				stagger = 3;
			if (end_pargroup == 0)
				p_radar.r_s.stagger = stagger;
			else
				p_slice[slice].r_s.stagger = stagger;
		}

		if (NULL != strstr(line, "<datatypez>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.CZ = 1;
			else
				p_radar.CZ = 0;
		}
		if (NULL != strstr(line, "<datatypeuz>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.UZ = 1;
			else
				p_radar.UZ = 0;
		}
		if (NULL != strstr(line, "<datatypev>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.V = 1;
			else
				p_radar.V = 0;
		}
		if (NULL != strstr(line, "<datatypew>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.W = 1;
			else
				p_radar.W = 0;
		}
		if (NULL != strstr(line, "<datatypezdr>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.ZDR = 1;
			else
				p_radar.ZDR = 0;
		}
		if (NULL != strstr(line, "<datatypephidp>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.PDP = 1;
			else
				p_radar.PDP = 0;
		}
		if (NULL != strstr(line, "<datatypekdp>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.KDP = 1;
			else
				p_radar.KDP = 0;
		}
		if (NULL != strstr(line, "<datatypeukdp>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.UKDP = 1;
			else
				p_radar.UKDP = 0;
		}
		if (NULL != strstr(line, "<datatyperhohv>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.RHV = 1;
			else
				p_radar.RHV = 0;
		}
		if (NULL != strstr(line, "<datatypeldr>"))
		{
			if (NULL != strstr(line, "On"))
				p_radar.LDR = 1;
			else
				p_radar.LDR = 0;
		}
		if (NULL != strstr(line, "<anglesync>"))
		{
			if (NULL != strstr(line, "On"))
				p_slice[slice].anglesync = 1;
			else
				p_slice[slice].anglesync = 0;
		}

		if (NULL != strstr(line, "<dyn"))
		{
			fill_dyn_info(line, &p_slice[slice]);
		}

		if (NULL != strstr(line, "<noise_power_dbz>"))
		{
			return_value(line, tmp_str);
			p_slice[slice].noiseh = (float) atof(tmp_str);
		}

		if (NULL != strstr(line, "<noise_power_dbz_dpv>"))
		{
			return_value(line, tmp_str);
			p_slice[slice].noisev = (float) atof(tmp_str);
		}
		if (NULL != strstr(line, "<rspdphradconst>"))
		{
			return_value(line, tmp_str);

			i=0;
			tmp_str2=strtok(tmp_str, " ");
			p_slice[slice].radconsth[i]= (float) atof(tmp_str2);
			while (tmp_str2 != NULL)
			{
				// printf("radar constant h: %f \r\n", p_slice[slice].radconsth[i]);
				// printf("i=%d \r\n", i);
				i++;
				tmp_str2=strtok(NULL, " ");
				if (tmp_str2 != NULL) p_slice[slice].radconsth[i]= (float) atof(tmp_str2);
			}
		}
		if (NULL != strstr(line, "<rspdpvradconst>"))
		{
			return_value(line, tmp_str);

			i=0;
			tmp_str2=strtok(tmp_str, " ");
			p_slice[slice].radconstv[i]= (float) atof(tmp_str2);
			while (tmp_str2 != NULL)
			{
				// printf("radar constant v: %f \r\n", p_slice[slice].radconstv[i]);
				// printf("i=%d \r\n", i);
				i++;
				tmp_str2=strtok(NULL, " ");
				if (tmp_str2 != NULL) p_slice[slice].radconstv[i]= (float) atof(tmp_str2);
			}
		}
		if (NULL != strstr(line, "<pw_index>"))
		{
			return_value(line, tmp_str);
			p_radar.pw_index = atoi(tmp_str);
			// printf("pulse width index: %d \r\n", p_radar.pw_index);
		}

		if (NULL != strstr(line, "<slicedata "))
		{
			strcpy(tmp_str, " time");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_slice[slice].time, tmp_str, SLICE_TIME_LEN);

			strcpy(tmp_str, "\" date");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_slice[slice].date, tmp_str, SLICE_DATE_LEN);
		}



		if (NULL != strstr(line, "<rayinfo "))
		{
			strcpy(tmp_str, "rays");
			return_info(line, tmp_str, tmp_str);
			p_slice[slice].rays = atoi(tmp_str);

			strcpy(tmp_str, "blobid");
			return_info(line, tmp_str, tmp_str);
			blobid = atoi(tmp_str);

			strcpy(tmp_str, "depth");
			return_info(line, tmp_str, tmp_str);
			blob_info[blobid].depth = atoi(tmp_str);

			strcpy(tmp_str, "refid");
			return_info(line, tmp_str, tmp_str);
			blob_info[blobid].type = BLOB_UNUSED;
			if (strcmp(tmp_str, "startangle") == 0) {
				blob_info[blobid].type = BLOB_STARTANGLE;
				p_slice[slice].blobid_startangle = blobid;
			}
			if (strcmp(tmp_str, "stopangle") == 0) {
				blob_info[blobid].type = BLOB_STOPANGLE;
				p_slice[slice].blobid_stopangle = blobid;
			}
			if (strcmp(tmp_str, "timestamp") == 0)
				blob_info[blobid].type = BLOB_TIMESTAMP;

			blob_info[blobid].slice = slice;
		}

		if (NULL != strstr(line, "<rawdata "))
		{
			strcpy(tmp_str, "rays");
			return_info(line, tmp_str, tmp_str);
			p_slice[slice].rays = atoi(tmp_str);

			strcpy(tmp_str, "bins");
			return_info(line, tmp_str, tmp_str);
			p_slice[slice].bins = atoi(tmp_str);

			strcpy(tmp_str, "blobid");
			return_info(line, tmp_str, tmp_str);
			blobid = atoi(tmp_str);
			p_slice[slice].blobid = blobid;

			strcpy(tmp_str, "depth");
			return_info(line, tmp_str, tmp_str);
			blob_info[blobid].depth = atoi(tmp_str);

			strcpy(tmp_str, "type");
			return_info(line, tmp_str, tmp_str);
				strncpy(p_slice[slice].typestr, tmp_str, SLICE_MAX_DATATYPE_LEN);

			blob_info[blobid].type = BLOB_UNUSED;
			if ( (strcmp(tmp_str, "dBZ") == 0) ||
				(strcmp(tmp_str, "dBZv") == 0))
				blob_info[blobid].type = BLOB_RAWDATACZ;
			if ( (strcmp(tmp_str, "dBuZ") == 0) ||
				(strcmp(tmp_str, "dBuZv") == 0))
				blob_info[blobid].type = BLOB_RAWDATAUZ;
			if ( (strcmp(tmp_str, "V") == 0) ||
				(strcmp(tmp_str, "Vv") == 0) ||
				(strcmp(tmp_str, "Vu") == 0) ||
				(strcmp(tmp_str, "Vvu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAV;
			if ( (strcmp(tmp_str, "W") == 0) ||
				(strcmp(tmp_str, "Wv") == 0) ||
				(strcmp(tmp_str, "Wu") == 0) ||
				(strcmp(tmp_str, "Wvu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAW;
			if ( (strcmp(tmp_str, "ZDR") == 0) ||
				(strcmp(tmp_str, "ZDRu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAZDR;
			if ( (strcmp(tmp_str, "PhiDP") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAPDP;
			if ( (strcmp(tmp_str, "uPhiDP") == 0) ||
				(strcmp(tmp_str, "uPhiDPu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAUPDP;
			if ( (strcmp(tmp_str, "KDP") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAKDP;
			if ( (strcmp(tmp_str, "uKDP") == 0) ||
				(strcmp(tmp_str, "uKDPu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATAUKDP;
			if ( (strcmp(tmp_str, "RhoHV") == 0) ||
				(strcmp(tmp_str, "RhoHVu") == 0) )
				blob_info[blobid].type = BLOB_RAWDATARHV;
			if (strcmp(tmp_str, "LDR") == 0)
			blob_info[blobid].type = BLOB_RAWDATALDR;
			if ( (strcmp(tmp_str, "SQI") == 0) ||
				(strcmp(tmp_str, "SQIv") == 0) ||
				(strcmp(tmp_str, "SQIu") == 0) ||
				(strcmp(tmp_str, "SQIvu") == 0) )
			blob_info[blobid].type = BLOB_RAWDATASQI;
			if ( (strcmp(tmp_str, "ISO0") == 0) )
				blob_info[blobid].type = BLOB_COSMOISO0;
			if ( (strcmp(tmp_str, "TEMP") == 0) )
				blob_info[blobid].type = BLOB_COSMOTEMP;
			if ( (strcmp(tmp_str, "COORD") == 0) )
				blob_info[blobid].type = BLOB_COORD;
			if ( (strcmp(tmp_str, "VIS") == 0) )
				blob_info[blobid].type = BLOB_VISIBILITY;

			blob_info[blobid].slice = slice;
		}

		if (NULL != strstr(line, "<sensorinfo ")) {
			strcpy(tmp_str, "type");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_sensor.type, tmp_str, MAX_SENSORTYPE_LEN);

			strcpy(tmp_str, "id");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_sensor.id, tmp_str, MAX_SENSORID_LEN);

			strcpy(tmp_str, "name");
			return_info(line, tmp_str, tmp_str);
			strncpy(p_sensor.name, tmp_str, MAX_SENSORNAME_LEN);
		}

		if (NULL != strstr(line, "<lat>"))
		{
			return_value(line, tmp_str);
			p_sensor.lat = atof(tmp_str);
		}
		if (NULL != strstr(line, "<lon>"))
		{
			return_value(line, tmp_str);
			p_sensor.lon = atof(tmp_str);
		}
		if (NULL != strstr(line, "<alt>"))
		{
			return_value(line, tmp_str);
			p_sensor.alt = atof(tmp_str);
		}
		if (NULL != strstr(line, "<wavelen>"))
		{
			return_value(line, tmp_str);
			p_sensor.wavelen = atof(tmp_str);
		}
		if (NULL != strstr(line, "<beamwidth>"))
		{
			return_value(line, tmp_str);
			p_sensor.beamwidth = atof(tmp_str);
		}
	}
	/* end reading file header  */

	if (p_radar.numele != (slice + 1)) {
		printf("Invalid number of slices (is %d, should be %d)\r\n",
			slice + 1, p_radar.numele);
		fclose(fp);
		return -1;
	}

	fgetpos(fp, &fileposition);
	strncpy(filename, fname, MAXFILENAME);
	fileheader_read = 1;

	fclose(fp);
	return 0;
}

static void fill_dyn_info(char *line, struct params_slice *slice)
{
  /*retrieve dyn information*/
  char type = 0;
  int i = 0, j = 0, num = 0;
  char read_number = 0;
  char number[2][8];

  memset(number, 0, 2 * 8);

  if (NULL != strstr(line, "dynz ")) type = 1;
  if (NULL != strstr(line, "dynv ")) type = 2;
  if (NULL != strstr(line, "dynw ")) type = 3;
  if (NULL != strstr(line, "dynzdr ")) type = 4;
  if (NULL != strstr(line, "dynldr ")) type = 5;
  if (NULL != strstr(line, "dynkdp ")) type = 6;


  /*there are always 2 numbers (min and max)*/
  while (line[i] != '>')
    {
      if (line[i] == '\"')
	{
	  j = 0;
	  if (read_number)
            {
	      /*1 number already read, add 1 to num*/
	      read_number = 0;
	      num++;
            }
	  else
            {
	      read_number = 1;
            }
	}
      else
	{
	  if (read_number)
            {
	      number[num][j++] = line[i];
            }
	}
      i++;
    }

  switch(type)
    {
    case 1:
      {
	slice->dynz[0] = (float) atof(number[0]);
	slice->dynz[1] = (float) atof(number[1]);
	break;
      }
    case 2:
      {
	slice->dynv[0] = (float) atof(number[0]);
	slice->dynv[1] = (float) atof(number[1]);
	break;
      }
    case 3:
      {
	slice->dynw[0] = (float) atof(number[0]);
	slice->dynw[1] = (float) atof(number[1]);
	break;
      }
    case 4:
      {
	slice->dynzdr[0] = (float) atof(number[0]);
	slice->dynzdr[1] = (float) atof(number[1]);
	break;
      }
    case 5:
      {
	slice->dynldr[0] = (float) atof(number[0]);
	slice->dynldr[1] = (float) atof(number[1]);
	break;
      }
    case 6:
      {
	slice->dynkdp[0] = (float) atof(number[0]);
	slice->dynkdp[1] = (float) atof(number[1]);
	break;
      }
    default:
      {
	break;
      }
    }
  return;
}

static void return_info(char *line, char *pattern, char *result)
{
  /*return values between " "*/
  char *tmp = NULL;
  int size_str = 0;
  int i = 0, j = 0;

  tmp = strstr(line, pattern);

  if (NULL == tmp)
    {
      return;
    }

  size_str = strlen(pattern);
  /*info starts at position strlen + 2 (to account for ' =" ') */
  i = size_str + 2;
  j = 0;
  memset(result, 0, MAX_LINE_SIZE);

  while (tmp[i] != '\"')
    result[j++] = tmp[i++];

  return;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
