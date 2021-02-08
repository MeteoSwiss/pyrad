/*
  *************************************************************

  IDL to C interface for rainbow files

  *************************************************************

  Filename:         idl_rainbow.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-19
  Last update:      2015-01-16 (fvj): added function to compress data
    
  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Mapping of idl external calls to the corresponding C functions.

  *************************************************************
  =============================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "idl_export.h"
#include "rainbow.h"

#ifdef __cplusplus
   extern "C" {
#endif

/* FUNCTIONS DEALING WITH DATA COMPRESSION */

/**
 * get expected size of compressed data
 *
 * @param nelem			number of elements in data to be compressed
 * @param datatype		idl data type
 * @param nbytescdata	expected size of compressed data 
 *
 * @return 			0 if successful or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getCompressedDataSize(int argc, void *argv[])
{
  IDL_LONG idl_type;
  IDL_LONG nelem;  
  unsigned long *nbytescdata;
  unsigned long nbytesudata;
  int ok;

  if (argc != 3) 
  {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  
  nelem    = *((IDL_LONG *)argv[0]); 
  idl_type = *((IDL_LONG *)argv[1]) & 0xffff;
  nbytescdata = ((IDL_ULONG *)argv[2]);
  
  switch (idl_type)
  {
    case IDL_TYP_BYTE:
	  nbytesudata=nelem;
	  break;
	case IDL_TYP_UINT:
	  nbytesudata=2*nelem;
	  break;
	case IDL_TYP_INT:
	  nbytesudata=2*nelem;
	  break;
	default:
	  printf("%s:%s: ERROR: Wrong datatype: type number %ld\r\n", 
	   __FILE__, __func__, idl_type);
	  return -1;
  }
  
  ok=rainbow_getCompressedDataSize(nbytesudata, nbytescdata);  
  if (ok <0)
  {
    printf("%s:%s: ERROR: Unable to determine expected compressed data size\r\n", __FILE__, __func__);
    return -1;
  }
  //printf("%s:%s: Expected compressed data size (bytes): %lu\r\n", __FILE__, __func__, *nbytescdata);
  return 0;  
}


/**
 * compress data using zlib library
 *
 * @param nelem		number of elements in data to be compressed
 * @param idl_type	idl data type
 * @param udata		pointer to uncompressed data
 * @param cdata		pointer where to place the compressed data 
 *
 * @return 			0 if successful or -1 if an error occurs -2 if assigned memory for array (in IDL program) too small
 */
IDL_LONG idl_rainbow_compressData(int argc, void *argv[])
{
  IDL_LONG idl_type;
  IDL_LONG nelem;
  unsigned char *udata;
  unsigned char *cdata;  
  unsigned long nbytesudata;
  unsigned long *nbytescdata;
  
  if (argc != 5) 
  {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  
  nelem    = *((IDL_LONG *)argv[0]); 
  idl_type = *((IDL_LONG *)argv[1]) & 0xffff;
  udata = (unsigned char *)argv[2];
  cdata = (unsigned char *)argv[3];
  nbytescdata = (IDL_ULONG *)argv[4];
  
  switch (idl_type)
  {
    case IDL_TYP_BYTE:
	  nbytesudata=nelem;
	  break;
	case IDL_TYP_UINT:
	  nbytesudata=2*nelem;
	  break;
	case IDL_TYP_INT:
	  nbytesudata=2*nelem;
	  break;
	default:
	  printf("%s:%s: ERROR: Wrong datatype: type number %ld\r\n", 
	   __FILE__, __func__, idl_type);
	  return -1;
  }

  if (!udata) {
    printf("%s:%s: ERROR: malloc()\r\n", __FILE__, __func__);
    return -1;
  }    
  return rainbow_compressData(nbytesudata, udata, nbytescdata, cdata);  
}

/* END OF FUNCTIONS DEALING WITH DATA COMPRESSION */


/* FUNCTIONS DEALING WITH SENSOR INFORMATION */

/**
 * Get sensor longitude wgs84
 *
 * @param filename Name of rainbow raw file
 
 * @return sensor longitude [Deg], -99999. otherwise
 */
float idl_rainbow_getLongitude(int argc, void *argv[])
{
  IDL_STRING fname;
  float lon;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  lon = rainbow_getLongitude(fname.s);
  if (lon == -99999.) {
    printf("%s:%s: Could not get sensor longitude file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -99999.;
  }  
  return lon;
}


/**
 * Get sensor latitude wgs84
 *
 * @param filename Name of rainbow raw file
 
 * @return sensor latitude [Deg], -99999. otherwise
 */
float idl_rainbow_getLatitude(int argc, void *argv[])
{
  IDL_STRING fname;
  float lat;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  lat = rainbow_getLatitude(fname.s);
  if (lat == -99999.) {
    printf("%s:%s: Could not get sensor latitude file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -99999.;
  }  
  return lat;
}

/**
 * Get sensor antenna height wgs84 [m msl]
 *
 * @param filename Name of rainbow raw file
 
 * @return sensor antenna height [m], -99999. otherwise
 */
float idl_rainbow_getAltitude(int argc, void *argv[])
{
  IDL_STRING fname;
  float alt;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  alt = rainbow_getAltitude(fname.s);
  if (alt == -99999.) {
    printf("%s:%s: Could not get sensor antenna height file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -99999.;
  }  
  return alt;
}

/**
 * Get sensor wavelength
 *
 * @param filename Name of rainbow raw file
 
 * @return sensor wavelength [m], -99999. otherwise
 */
float idl_rainbow_getWavelength(int argc, void *argv[])
{
  IDL_STRING fname;
  float wavelen;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  wavelen = rainbow_getWavelength(fname.s);
  if (wavelen == -99999.) {
    printf("%s:%s: Could not get sensor wavelength file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -99999.;
  }  
  return wavelen;
}

/**
 * Get sensor beamwidth
 *
 * @param filename Name of rainbow raw file
 
 * @return sensor beam width [deg], -99999. otherwise
 */
float idl_rainbow_getBeamwidth(int argc, void *argv[])
{
  IDL_STRING fname;
  float beamwidth;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  beamwidth = rainbow_getBeamwidth(fname.s);
  if (beamwidth == -99999.) {
    printf("%s:%s: Could not get sensor beam width file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -99999.;
  }  
  return beamwidth;
}

/**
 * Get the sensor type from a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file  
 * @param str      Empty string to copy the sensor type
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSensorType(int argc, void *argv[])
{  
  IDL_STRING fname;
  IDL_STRING type_str;
  char *str = NULL;
  int ok;
  unsigned long *slen;

  if (argc != 3) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);  
  type_str = *((IDL_STRING *)argv[1]);
  slen = (IDL_ULONG *)argv[2];
  

  ok = rainbow_getSensorType(fname.s, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the sensor type from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  //printf("%s:%s: Sensor type %s\r\n", __FILE__, __func__, str);
  strncpy(type_str.s, str, type_str.slen);
  *slen=(unsigned long) strlen(str);  
  //printf("%s:%s: Sensor type string length %lu\r\n", __FILE__, __func__, *slen);
  return 0;
}

/**
 * Get the sensor ID from a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file  
 * @param str      Empty string to copy the sensor ID
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSensorID(int argc, void *argv[])
{  
  IDL_STRING fname;
  IDL_STRING id_str;
  char *str = NULL;
  int ok;
  unsigned long *slen;

  if (argc != 3) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);  
  id_str = *((IDL_STRING *)argv[1]);
  slen = (IDL_ULONG *)argv[2];

  ok = rainbow_getSensorID(fname.s, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the sensor ID from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  //printf("%s:%s: Sensor id %s\r\n", __FILE__, __func__, str);
  strncpy(id_str.s, str, id_str.slen);
  *slen=(unsigned long) strlen(str);  
  return 0;
}

/**
 * Get the sensor name from a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file  
 * @param str      Empty string to copy the sensor name
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSensorName(int argc, void *argv[])
{  
  IDL_STRING fname;
  IDL_STRING name_str;
  char *str = NULL;
  int ok;
  unsigned long *slen;

  if (argc != 3) {
    printf("%s:%s: ERROR: Wrong number of arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);  
  name_str = *((IDL_STRING *)argv[1]);
  slen = (IDL_ULONG *)argv[2];

  ok = rainbow_getSensorName(fname.s, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the sensor name from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  //printf("%s:%s: Sensor name %s\r\n", __FILE__, __func__, str);
  strncpy(name_str.s, str, name_str.slen);
  *slen=(unsigned long) strlen(str);  
  return 0;
}

/* END OF FUNCTIONS DEALING WITH SENSOR INFORMATION */


/**
 * Get the horizontal channel noise power [dBZ at 1 km]
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 *
 * @return noise power [dBZ at 1 km], -99999, otherwise
 */
float idl_rainbow_getNoisePowerZh(int argc, void *argv[])
{
  IDL_STRING fname;
  IDL_LONG slicenum;
  float NoisePowerZh;
  
  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  
  NoisePowerZh=rainbow_getNoisePowerZh(fname.s, (int)slicenum);
  
  if (NoisePowerZh == -99999.) {
    printf("%s:%s: Could not get the horizontal channel noise power of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -99999.;
  }
  return NoisePowerZh;
}

/**
 * Get the vertical channel noise power [dBZ at 1 km]
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 *
 * @return noise power [dBZ at 1 km], -99999, otherwise
 */
float idl_rainbow_getNoisePowerZv(int argc, void *argv[])
{
  IDL_STRING fname;
  IDL_LONG slicenum;
  float NoisePowerZv;
  
  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -99999.;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  
  NoisePowerZv=rainbow_getNoisePowerZv(fname.s, (int)slicenum);
  
  if (NoisePowerZv == -99999.) {
    printf("%s:%s: Could not get the horizontal channel noise power of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -99999.;
  }
  return NoisePowerZv;
}

/**
 * Get the horizontal channel radar constants
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param radarconstant pointer to radar constants array (float array) to write the radar constants
 * @param radarconstant_len length of radar cosntants array
 *
 * @return 0 on success, -1 otherwise
 */
float idl_rainbow_getRadarConstanth(int argc, void *argv[])
{
  IDL_STRING fname;
  IDL_LONG slicenum;
  float *radarconstant;
  unsigned long radarconstant_len;
  int ok;  
  
  if (argc != 4) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  radarconstant = (float *)argv[2];
  radarconstant_len = *(IDL_LONG *)argv[3];
  
  ok = rainbow_getRadarConstanth(fname.s, (int)slicenum, radarconstant, radarconstant_len);
  if (ok < 0) {
    printf("%s:%s: Could not get the H channel radar constants from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }    
  return 0;
}

/**
 * Get the horizontal channel radar constants
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param radarconstant pointer to radar constants array (float array) to write the radar constants
 * @param radarconstant_len length of radar cosntants array
 *
 * @return 0 on success, -1 otherwise
 */
float idl_rainbow_getRadarConstantv(int argc, void *argv[])
{
  IDL_STRING fname;
  IDL_LONG slicenum;
  float *radarconstant;
  unsigned long radarconstant_len;
  int ok;  
  
  if (argc != 4) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  radarconstant = (float *)argv[2];
  radarconstant_len = *(IDL_LONG *)argv[3];
  
  ok = rainbow_getRadarConstantv(fname.s, (int)slicenum, radarconstant, radarconstant_len);
  if (ok < 0) {
    printf("%s:%s: Could not get the V channel radar constants from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }    
  return 0;
}

/**
 * Get pulse width index stored in a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 *
 * @return pulse width index, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getPulseWidthIndex(int argc, void *argv[])
{
  int ind_pw;
  IDL_STRING fname;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  ind_pw = rainbow_getPulseWidthIndex(fname.s);
  if (ind_pw < 0) {
    printf("%s:%s: Could not get the pulsewidth index from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }

  return (IDL_LONG)ind_pw;
}






/**
 * Get number of slices stored a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 *
 * @return number of slices, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getNumberOfSlices(int argc, void *argv[])
{
  int numslices;
  IDL_STRING fname;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  numslices = rainbow_getNumberOfSlices(fname.s);
  if (numslices < 0) {
    printf("%s:%s: Could not get the number of slices from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }

  //printf("Number of slices: %d\r\n", numslices);

  return (IDL_LONG)numslices;
}

/**
 * Get the scan name of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param str      Empty string to copy the scan name
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getScanName(int argc, void *argv[])
{
  IDL_STRING fname;
  IDL_STRING scanname_str;
  char *str = NULL;
  int ok;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  scanname_str = *((IDL_STRING *)argv[1]);

  ok = rainbow_getScanName(fname.s, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the scan name from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  strncpy(scanname_str.s, str, scanname_str.slen);
  return 0;
}

/**
 * Get the range resolution of a scan (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 *
 * @return range resolution in m, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getRangeResolution(int argc, void *argv[])
{
  float rangeres;
  IDL_STRING fname;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  rangeres = rainbow_getRangeResolution(fname.s); /* [km] */
  if (rangeres < 0) {
    printf("%s:%s: Could not get range resolution file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  return (IDL_LONG)(rangeres * 1000); /* convert [km] to [m] */
}

/**
 * Get the range sampling of a scan (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 *
 * @return range sampling, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getRangeSampling(int argc, void *argv[])
{
  int rangesamp;
  IDL_STRING fname;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  rangesamp = rainbow_getRangeSampling(fname.s);
  if (rangesamp < 0) {
    printf("%s:%s: Could not get range sampling file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  return (IDL_LONG)(rangesamp);
}

/**
 * Get a list of slice angles (= elevation angles for a volumne scan)
 *
 * @param filename Name of rainbow raw file
 * @param sliceangles pointer to slice angle array (float array) to write the slice angles
 * @param sliceangles_len length of slice array
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSliceAngles(int argc, void *argv[])
{
  IDL_STRING fname;
  float *sliceangles;
  IDL_LONG sliceangles_len;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  sliceangles = (float *)argv[1];
  sliceangles_len = *(IDL_LONG *)argv[2];

  if (rainbow_getSliceAngles(fname.s, sliceangles, sliceangles_len) < 0) {
    printf("%s:%s: Could not get slice angles from file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1;
  }
  return 0;
}




/**
 * Get angle step
 *
 * @param filename Name of rainbow raw file
 
 * @return angle step [Deg], -1 otherwise
 */
float idl_rainbow_getAngleStep(int argc, void *argv[])
{
  IDL_STRING fname;
  float anglestep;
  
  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  anglestep = rainbow_getAngleStep(fname.s);
  if (anglestep < 0.) {
    printf("%s:%s: Could not get angle step file \'%s\'\r\n",
	   __FILE__, __func__, fname.s);
    return -1.;
  }
  //printf("%s:%s: Angle step: %f\r\n", __FILE__, __func__, anglestep);
  return anglestep;
}


/**
 * Get the fixed angle of a slice of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 * @param pointer to angle variable to write (float)
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getFixedAngle(int argc, void *argv[])
{
  IDL_LONG slicenum;
  IDL_STRING fname;
  float *fixedangle;
  int ok;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  fixedangle = (float *)argv[2];
  
  ok = rainbow_getFixedAngle(fname.s, (int)slicenum, fixedangle);
  if (ok < 0) {
    printf("%s:%s: Could not get the fixedangle of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }
  return 0;
}

/**
 * Get the time of a slices of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 * @param str      Empty string to copy the time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSliceTime(int argc, void *argv[])
{
  IDL_LONG slicenum;
  IDL_STRING fname;
  IDL_STRING time_str;
  char *str = NULL;
  int ok;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  time_str = *((IDL_STRING *)argv[2]);

  ok = rainbow_getSliceTime(fname.s, (int)slicenum, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the acquisition time of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }
  strncpy(time_str.s, str, time_str.slen);
  return 0;
}

/**
 * Get the antenna speed during the slice recording
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 * @param str      Empty string to copy the time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getAntSpeed(int argc, void *argv[])
{
  IDL_LONG slicenum;
  IDL_STRING fname;
  float *antspeed;
  int ok;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  antspeed = (float *)argv[2];

  ok = rainbow_getAntSpeed(fname.s, (int)slicenum, antspeed);
  if (ok < 0) {
    printf("%s:%s: Could not get the antenna speed of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }
  return 0;
}

/**
 * Get the date of a slices of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 * @param str      Empty string to copy the time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSliceDate(int argc, void *argv[])
{
  IDL_LONG slicenum;
  IDL_STRING fname;
  IDL_STRING date_str;
  char *str = NULL;
  int ok;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  date_str = *((IDL_STRING *)argv[2]);

  ok = rainbow_getSliceDate(fname.s, (int)slicenum, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the acquisition date of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }
  strncpy(date_str.s, str, date_str.slen);
  return 0;
}

/**
 * Get the datatype a slices of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 * @param str      Empty string to copy the time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_rainbow_getSliceDatatype(int argc, void *argv[])
{
  IDL_LONG slicenum;
  IDL_STRING fname;
  IDL_STRING type_str;
  char *str = NULL;
  int ok;
  int len;
  int k;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  type_str = *((IDL_STRING *)argv[2]);

  ok = rainbow_getSliceDatatype(fname.s, (int)slicenum, &str);
  if (ok < 0) {
    printf("%s:%s: Could not get the acquisition date of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }
  strncpy(type_str.s, str, type_str.slen);
  len = strlen(str);
  for (k=len;k<type_str.slen;k++) {
    type_str.s[k] = ' ';
  }
  return 0;
}

/**
 * Get number of angles of a slices of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 *
 * @return number of angles, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getNumberOfAngles(int argc, void *argv[])
{
  int numangles;
  IDL_LONG slicenum;
  IDL_STRING fname;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);

  numangles = rainbow_getNumberOfAngles(fname.s, (int)slicenum);
  if (numangles < 0) {
    printf("%s:%s: Could not get the number of angles of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }

  //printf("Number of angles: %d\r\n", numangles);

  return (IDL_LONG)numangles;
}

/**
 * Get number of range bins of a slices of a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file 
 * @param slicenum Number of the slice (0=first slice)
 *
 * @return number of bins, or -1 if an error occurs
 */
IDL_LONG idl_rainbow_getNumberOfRangeBins(int argc, void *argv[])
{
  int numbins;
  IDL_LONG slicenum;
  IDL_STRING fname;

  if (argc != 2) {
    printf(__FILE__ " ERROR: Not enough arguments\r\n");
    return -1;
  }

  fname = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);

  numbins = rainbow_getNumberOfRangeBins(fname.s, (int)slicenum);
  if (numbins < 0) {
    printf("%s:%s: Could not get the number of angles of slice %ld from file \'%s\'\r\n",
	   __FILE__, __func__, slicenum, fname.s);
    return -1;
  }

  //printf("Number of bins: %d\r\n", numbins);

  return (IDL_LONG)numbins;
}

/**
 * Get data of a slice stored a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param ncols
 * @param nrows
 * @param idl_type
 * @param data
 *
 * @return 0 on sucess, -1 otherwise
 */
IDL_LONG idl_rainbow_getSliceData(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG slicenum;
  IDL_LONG ncols;
  IDL_LONG nrows;
  IDL_LONG idl_type;
  float *data;
  int ret;

  if (argc != 6) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  ncols = *((IDL_LONG *)argv[2]); // angles
  nrows = *((IDL_LONG *)argv[3]); // bins
  idl_type = *((IDL_LONG *)argv[4]) & 0xffff;
  data = ((float *)argv[5]);

  if (idl_type != IDL_TYP_FLOAT) {
    printf("%s:%s: ERROR: Wrong datatype: type number %ld instead of %d\r\n", 
	   __FILE__, __func__, idl_type, IDL_TYP_LONG);
    return -1;
  }

  if (!data) {
    printf("%s:%s: ERROR: malloc()\r\n", __FILE__, __func__);
    return -1;
  }

  //printf("Get data from slice %d\r\n", slicenum);

  ret = rainbow_getSliceData(filename.s, slicenum, data, (long)(ncols*nrows));
  return (IDL_LONG)ret;
}


/**
 * Get start angles of a slice stored in a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param nangs    Number of angles
 * @param idl_type Type of array
 * @param data     Pointer to data
 *
 * @return 0 on sucess, -1 otherwise
 */
IDL_LONG idl_rainbow_getStartAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG slicenum;
  IDL_LONG nangs;
  IDL_LONG idl_type;
  float *data;
  int ret;

  if (argc != 5) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  nangs    = *((IDL_LONG *)argv[2]); // number of angles
  idl_type = *((IDL_LONG *)argv[3]) & 0xffff;
  data = ((float *)argv[4]);

  if (idl_type != IDL_TYP_FLOAT) {
    printf("%s:%s: ERROR: Wrong datatype: type number %ld instead of %d\r\n", 
	   __FILE__, __func__, idl_type, IDL_TYP_LONG);
    return -1;
  }

  if (!data) {
    printf("%s:%s: ERROR: malloc()\r\n", __FILE__, __func__);
    return -1;
  }

  ret = rainbow_getStartAngles(filename.s, slicenum, data, (long)nangs);
  return (IDL_LONG)ret;
}

/**
 * Get start angles of a slice stored in a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param nangs    Number of angles
 * @param idl_type Type of array
 * @param data     Pointer to data
 *
 * @return 0 on sucess, -1 otherwise
 */
IDL_LONG idl_rainbow_getAngles(int argc, void *argv[])
{
  return idl_rainbow_getStartAngles(argc, argv);
}

/**
 * Get stop angles of a slice stored in a rainbow raw file (either .vol, 
 * .azi or .ele file)
 *
 * @param filename Name of rainbow raw file
 * @param slicenum Number of the slice (0=first slice)
 * @param nangs    Number of angles
 * @param idl_type Type of array
 * @param data     Pointer to data
 *
 * @return 0 on sucess, -1 otherwise
 */
IDL_LONG idl_rainbow_getStopAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG slicenum;
  IDL_LONG nangs;
  IDL_LONG idl_type;
  float *data;
  int ret;

  if (argc != 5) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  slicenum = *((IDL_LONG *)argv[1]);
  nangs    = *((IDL_LONG *)argv[2]); // number of angles
  idl_type = *((IDL_LONG *)argv[3]) & 0xffff;
  data = ((float *)argv[4]);

  if (idl_type != IDL_TYP_FLOAT) {
    printf("%s:%s: ERROR: Wrong datatype: type number %ld instead of %d\r\n", 
	   __FILE__, __func__, idl_type, IDL_TYP_LONG);
    return -1;
  }

  if (!data) {
    printf("%s:%s: ERROR: malloc()\r\n", __FILE__, __func__);
    return -1;
  }

  ret = rainbow_getStopAngles(filename.s, slicenum, data, (long)nangs);
  return (IDL_LONG)ret;
}

#ifdef __cplusplus
   } /* extern "C" */
#endif
