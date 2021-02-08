/*
  *************************************************************

  Rainbow 

  *************************************************************

  Filename:         rainbow.h 
  Author:           Andreas Leuenberger
  Creation date:    2012-10-19
  Last update:      2013-03-08
    
  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------

  *************************************************************
  =============================================================
 */

#ifndef __RAINBOW_H
#define __RAINBOW_H

#ifdef __cplusplus
extern "C" {
#endif

/* General file info */
int rainbow_getNumberOfSlices(char *fname);
float rainbow_getRangeResolution(char *fname);
int rainbow_getRangeSampling(char *fname);

int rainbow_getSliceAngles(char *fname, float *sliceangles, int len);
float rainbow_getAngleStep(char *fname);

/* Sensor info */
float rainbow_getLongitude(char *fname);
float rainbow_getLatitude(char *fname);
float rainbow_getAltitude(char *fname);
float rainbow_getWavelength(char *fname);
float rainbow_getBeamwidth(char *fname);
int rainbow_getSensorType(char *fname, char **str);
int rainbow_getSensorID(char *fname, char **str);
int rainbow_getSensorName(char *fname, char **str);

/* Slice info */
int rainbow_getNumberOfAngles(char *fname, int slicenum);
int rainbow_getScanName(char *fname, char **str);
int rainbow_getNumberOfRangeBins(char *fname, int slicenum);
int rainbow_getSliceDate(char *fname, int slicenum, char **str);
int rainbow_getSliceTime(char *fname, int slicenum, char **str);
int rainbow_getSliceDatatype(char *fname, int slicenum, char **str);
int rainbow_getFixedAngle(char *fname, int slicenum, float *elevation);
int rainbow_getStartAngles(char *fname, int slicenum, float *data, long nangles);
int rainbow_getStopAngles(char *fname, int slicenum, float *data, long nangles);
int rainbow_getSliceData(char *fname, int slicenum, float *data, long nvals);
int rainbow_getAntSpeed(char *fname, int slicenum, float *antspeed);

float rainbow_getNoisePowerZh(char *fname, int slicenum);
float rainbow_getNoisePowerZv(char *fname, int slicenum);

int rainbow_getRadarConstanth(char *fname, int slicenum, float *radarconstants, int len);
int rainbow_getRadarConstantv(char *fname, int slicenum, float *radarconstants, int len);
int rainbow_getPulseWidthIndex(char *fname);

/* Compression of parameters to put into BLOB */
int rainbow_compressData(const unsigned long nbytesudata, const unsigned char *udata, unsigned long *nbytescdata, unsigned char *cdata); // fvj 15.01.2015
int rainbow_getCompressedDataSize(const unsigned long nbytesudata, unsigned long *nbytescdata);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __RAINBOW_H */
