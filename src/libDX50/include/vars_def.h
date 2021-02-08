/*
  *************************************************************

  Definitions for rainbow file reading

  *************************************************************

  Filename:         vars_def.h 
  Author:           
  Creation date:    
  Last update:      2015-05-15 (fvj: added definitions for visibility data)
    
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

#ifndef __VARS_DEF_H
#define __VARS_DEF_H

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_LINE_SIZE 1024
#define MAX_SLICES 64
#define MAX_BLOBS 1024

#define BLOB_UNUSED      0
#define BLOB_RAWDATAUZ   1
#define BLOB_RAWDATACZ   2
#define BLOB_RAWDATAV    3
#define BLOB_RAWDATAW    4
#define BLOB_RAWDATAZDR  5
#define BLOB_RAWDATAPDP  6
#define BLOB_RAWDATAUPDP 7
#define BLOB_RAWDATAKDP  8
#define BLOB_RAWDATAUKDP 9
#define BLOB_RAWDATARHV  10
#define BLOB_RAWDATASQI  11
#define BLOB_RAWDATALDR  12
#define BLOB_STARTANGLE  13
#define BLOB_STOPANGLE   14
#define BLOB_TIMESTAMP   15
#define BLOB_NUMPULESE   16
#define BLOB_DATAFLAG    17
#define BLOB_TXPOWER     18
#define BLOB_COSMOISO0	 19
#define BLOB_COSMOTEMP	 20
#define BLOB_COORD	 	 21
#define BLOB_VISIBILITY	 22

#define SCAN_VOL 0
#define SCAN_PPI 1
#define SCAN_RHI 2

#define SLICE_TIME_LEN 9
#define SLICE_DATE_LEN 11

#define SLICE_MAX_DATATYPE_LEN 8

#define MAX_SCANNAME_LEN 128

#define MAX_SENSORTYPE_LEN 128
#define MAX_SENSORID_LEN 128
#define MAX_SENSORNAME_LEN 128

# define MAX_PULSE_WIDTH 8

struct blob_ray_info
   {
   /*bits per bin*/
   unsigned short int depth;
   /*slice id*/
   unsigned short int slice;
   /*type - refer to #defines above*/
   char type;
   char ruf[3];
   };

struct common_radar_slice
   {
   /*stop range in Km*/
   float stoprange;
   /*range res in Km*/
   float rangestep;
   /*range of first gate*/
   float startrange;
   /*PRFs (Hz)*/
   float highprf;
   float lowprf;
   /*antenna speed (deg/s)*/
   float antspeed;
   /*angle step (deg)*/
   float anglestep;
   /*samples in time*/
   float timesamp;
   /*number of samples per gate*/
   unsigned short int rangesamp;
   /*staggering ratio: stagger / (stagger - 1)*/
   unsigned short int stagger;
   };


struct params_radar
   {
   struct common_radar_slice r_s;
   /*number of elevations*/
   unsigned short int numele;
   char scanname[MAX_SCANNAME_LEN];
   /*data types used in scan*/
   char UZ:1;
   char CZ:1;
   char V:1;
   char W:1;
   char ZDR:1;
   char PDP:1;
   char KDP:1;
   char UKDP:1;
   char RHV:1;
   char LDR:1;
   char alignment_char:6;
   /* pulse width index */
   int pw_index;
   };



struct params_slice
   {
   struct common_radar_slice r_s;
   /*dynamic ranges - index 0 is min and 1 is max*/
   float dynz[2];
   float dynv[2];
   float dynw[2];
   float dynzdr[2];
   float dynldr[2];
   float dynkdp[2];
   /*elevation (vol or ppi scan) or azimuth (rhi) angle*/
   float angle;
   /*number of rays for this slice*/
   unsigned short int rays;
   /*number of bins for this slice*/
   unsigned short int bins;
   /*AngleSync*/
   /*if 1 samples in fixed angles (usually we use fixed angles)*/
   /*if 0 samples in fixed time*/
   char anglesync;
   /* data blobid of slice */
   short int blobid;
   /* startangles blobid of slice */
   short int blobid_startangle;
   /* stopangles blobid of slice */
   short int blobid_stopangle;
   /* date of slice 'YYYY-MM-DD' */
   char date[SLICE_DATE_LEN];
   /* time of slice 'HH:MM:SS' */
   char time[SLICE_TIME_LEN];
   /* data type of slice */
   char typestr[SLICE_MAX_DATATYPE_LEN];
   /* horizontal channel noise power [dBZ at 1 km] */
   float noiseh;
   /* vertical channel noise power [dBZ at 1 km] */
   float noisev;
   /* radar constants of the H channel, 1 per pulse width */
   float radconsth[MAX_PULSE_WIDTH];
   /* radar constants of the V channel, 1 per pulse width */
   float radconstv[MAX_PULSE_WIDTH];
   };

struct params_sensor
   {
   char type[MAX_SENSORTYPE_LEN];
   char id[MAX_SENSORID_LEN];
   char name[MAX_SENSORNAME_LEN];
   double lon;
   double lat;
   double alt;
   double wavelen;
   double beamwidth;
   };

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __VARS_DEF_H */
