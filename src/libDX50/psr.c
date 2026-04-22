/*
*************************************************************

 Power Spectrum Recording

*************************************************************

 Filename:         psr.c
 Author:           Andreas Leuenberger
 Creation date:    2012-10-25
 Last update:      2014-03-26

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

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "endianness.h"
#include "psr_header.h"
#include "psr.h"
#include "idl_export.h"

__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");

#ifdef __cplusplus
extern "C" {
#endif

  /* forward definition */
  static int readItemTable(struct psr_file *ptr, FILE *fh);
  static int getValueArrayFromCpiHeader(struct psr_file *ptr, int valueoffset, void *data);
  static int getValueFromCpiHeader(struct psr_file *ptr, int item, int valueoffset, void *data);
  static void *getCpiBase(struct psr_file *ptr, int itemnum);
  static void psr_file_clean(struct psr_file *ptr);
  static void printCpiHeader(void *pcpi);

  /* static variables */
  static struct psr_file psrfile;

  static void *decompress_buffer;
  static long decompress_len;
  static long decompress_position;

  /*
   * ========================================================================================
   * ----------------------- Version --------------------------------------------------------
   * ========================================================================================
   */

  /**
   * Print the version string to stdout.
   */
  void psr_printVersion(void)
  {
    printf("%s\r\n", PSR_VERSION_STR);
  }

  /*
   * ========================================================================================
   * ----------------------- PSR file header info -------------------------------------------
   * ========================================================================================
   */

  /**
   * Get header length
   *
   * @param fname  Name of psr file
   *
   * @return Length of psr file header in bytes, -1 on error.
   */
  long psr_getHeaderLength(char *fname)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return psrfile.header.headerlen;
  }

  /**
   * Get number of items from psr file
   *
   * @param fname  Name of psr file
   *
   * @return Number of psr items, or -1 on error.
   */
  long psr_getNumberOfItems(char *fname)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return psrfile.header.itemtypepsr;
  }

  /**
   * Get range start
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRangeStart(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.rspstart;
    return 0;
  }

  /**
   * Get range stop
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRangeStop(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.rspstop;
    return 0;
  }

  /**
   * Get range step
   *
   * @param fname  Name of psr file
   * @param value  Range step [km]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRangeStep(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.rspstep;
    return 0;
  }

  /**
   * Get two way path atmospheric attenuation [db/km]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTwoWayPathAttenuation(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.pathatt;
    return 0;
  }

  /**
   * Get start time of slice
   *
   * @param fname  Name of psr file
   * @param value  Time string
   * @param len    String len (without \0)
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getStartTime(char *fname, char *value, int len)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    strncpy(value, psrfile.header.datastarttime, len+1);
    if (len+1 > 0) {
      value[len] = '\n';
    }
    return 0;
  }

  /**
   * Get stop time of slice
   *
   * @param fname  Name of psr file
   * @param value  Time string
   * @param len    String len (without \0)
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getStopTime(char *fname, char *value, int len)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    strncpy(value, psrfile.header.datastoptime, len+1);
    if (len+1 > 0) {
      value[len] = '\n';
    }
    return 0;
  }

  /**
   * Get receiver loss of horizontal channel [db]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRxLossH(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.rxloss_h;
    return 0;
  }

  /**
   * Get receiver loss of vertical channel [db]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRxLossV(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.rxloss_v;
    return 0;
  }

  /**
   * Get transmitter loss of horizontal channel [db]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxLossH(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.txloss_h;
    return 0;
  }

  /**
   * Get transmitter loss of vertical channel [db]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxLossV(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.txloss_v;
    return 0;
  }

  /**
   * Get radome loss [db]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRadomeLoss(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.radomeloss;
    return 0;
  }

  /**
   * Get transmit power split factor of horizontal channel [1]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxPowerSplitFactorH(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.txpowerfact_h;
    return 0;
  }

  /**
   * Get transmit power split factor of vertical channel [1]
   *
   * @param fname  Name of psr file
   * @param value  Pointer to value to set
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxPowerSplitFactorV(char *fname, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    *value =  psrfile.header.txpowerfact_v;
    return 0;
  }

  /**
   * Get pulse width index
   *
   * @param fname  Name of psr file
   *
   * @return pulse width index
   *               0: short pulse 0.33 us
   *               1: short pulse 0.5 us
   *               2: short pulse 1.2 us
   *               3: short pulse 2.0 us
   *              -1: error
   */
  long psr_getPWIndex(char *fname)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    return psrfile.header.pwidth_id;
  }

  /**
   * Get horizontal noise power of a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getNoisePowerH(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.noisepwr_h[pw];
    return 0;
  }

  /**
   * Get vertical noise power of a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getNoisePowerV(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.noisepwr_v[pw];
    return 0;
  }

  /**
   * Get horizontal radar constant of a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRadarConstantH(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.radconst_h[pw];
    return 0;
  }

  /**
   * Get vertical radar constant of a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getRadarConstantV(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.radconst_v[pw];
    return 0;
  }

  /**
   * Get horizontal dbADU to dBm offset a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getdBadu2dBmOffsetH(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.adu2dbmoffset_h[pw];
    return 0;
  }

  /**
   * Get vertical dbADU to dBm offset a pulse width
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getdBadu2dBmOffsetV(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.adu2dbmoffset_v[pw];
    return 0;
  }

  /**
   * Get matched filter loss of a pulse width [db]
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getMatchedFilterLoss(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.mfloss[pw];
    return 0;
  }

  /**
   * Get matched filter bandwidth of a pulse width [kHz]
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getMatchedFilterBandwidth(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.mfbandwidth[pw];
    return 0;
  }

  /**
   * Get tx norm power [kW]
   *
   * @param fname  Name of psr file
   * @param pw     Pulse width index (0-3)
   * @param value  Pointer to return value
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxNormPower(char *fname, long pw, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if ((pw >= NUM_PW) || (pw < 0)) {
      return -1;
    }
    *value =  psrfile.header.txnormpower[pw];
    return 0;
  }

  /*
   * ========================================================================================
   * ----------------------- PSR header info of all items -----------------------------------
   * ========================================================================================
   */

  /**
   * Get number of items from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to az angles to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getAzimuthStartAngles(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_AZIMUTHSTART_OFFSET, (void *)data);
  }

  /**
   * Get number of items from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to az angles to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getAzimuthStopAngles(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_AZIMUTHSTOP_OFFSET, (void *)data);
  }

  /**
   * Get number of items from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to el angles to write
   * @return 0 on success or -1 on error.
   */
  long psr_getElevationStartAngles(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_ELEVATIONSTART_OFFSET, (void *)data);
  }

  /**
   * Get number of items from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to el angles to write
   * @return 0 on success or -1 on error.
   */
  long psr_getElevationStopAngles(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_ELEVATIONSTOP_OFFSET, (void *)data);
  }

  /**
   * Get an array with the number of tx pulses from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to data to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getNumberOfTxPulsesArray(char *fname, long *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_NUM_TXPULSES_OFFSET, (void *)data);
  }

  /**
   * Get an array with PRF frequencies from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to data to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getPrfArray(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_PRF_OFFSET, (void *)data);
  }

  /**
   * Get an array with the number of range gates from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to data to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getNumberOfRangeGatesArray(char *fname, long *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_NUM_RANGEGATES_OFFSET, (void *)data);
  }

  /**
   * Get an array with averaged tx pulse powerPRF from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to data to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getTxPowerArray(char *fname, float *data)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    return getValueArrayFromCpiHeader(&psrfile, CPI_AVGTXPOWER_OFFSET, (void *)data);
  }

  /**
   * Get an array with noise values from psr file
   *
   * @param fname  Name of psr file
   * @param data   Pointer to data to write
   *
   * @return 0 on success or -1 on error.
   */
  long psr_getNoiseArray(char *fname, float *data)
  {
    int k;

    if (psr_openFile(fname) < 0)
      return -1;

    if (psrfile.header.reduced) {
      k=0;
      while(k < psrfile.header.itemtypepsr) {
        data[k] = psrfile.item_table[k].noise;
        k++;
      }
    } else {
      fprintf(stderr, "%s:%s(): WARNING: No noise values in file\r\n",
              __FILE__, __func__);
      return -1;
    }
    return 0;
  }

  /*
   * ========================================================================================
   * ----------------------- PSR CPI header info --------------------------------------------
   * ========================================================================================
   */

  /**
   * Print CPI header
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   *
   * @return 0 on sucess, or -1 on error.
   */
  long psr_printCpiHeader(char *fname, int item)
  {
    void *pcpi;

    if (psr_openFile(fname) < 0)
      return -1;

    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }

    pcpi = getCpiBase(&psrfile, item);
    printf("CPI Header of Item %d\r\n", item);
    printCpiHeader(pcpi);

    return 0;
  }

  /**
   * Get number of tx pulses
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   *
   * @return Number of tx pulses of selected item, or -1 on error.
   */
  long psr_getNumberOfTxPulses(char *fname, int item)
  {
    long data = -1;

    if (psr_openFile(fname) < 0)
      return -1;

    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }

    getValueFromCpiHeader(&psrfile, item, CPI_NUM_TXPULSES_OFFSET, &data);

    return data;
  }

  /**
   * Get number of range gates
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   *
   * @return Number of range gates of selected item, or -1 on error.
   */
  long psr_getNumberOfRangeGates(char *fname, int item)
  {
    long data = -1;

    if (psr_openFile(fname) < 0)
      return -1;

    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }

    getValueFromCpiHeader(&psrfile, item, CPI_NUM_RANGEGATES_OFFSET, &data);

    return data;
  }

  /**
   * Get azimuth start angle
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  Azimuth start angle [deg]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getAzimuthStartAngle(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_AZIMUTHSTART_OFFSET, value);
    return 0;
  }

  /**
   * Get azimuth stop angle
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  Azimuth start angle [deg]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getAzimuthStopAngle(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_AZIMUTHSTOP_OFFSET, value);
    return 0;
  }

  /**
   * Get elevation start angle
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  Elevation start angle [deg]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getElevationStartAngle(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_ELEVATIONSTART_OFFSET, value);
    return 0;
  }

  /**
   * Get elevation stop angle
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  Elevation start angle [deg]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getElevationStopAngle(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_ELEVATIONSTOP_OFFSET, value);
    return 0;
  }

  /**
   * Get prf
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  prf [Hz]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getPrf(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_PRF_OFFSET, value);
    return 0;
  }

  /**
   * Get averaged TX power
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  tx power [kW]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getTxPower(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_AVGTXPOWER_OFFSET, value);
    return 0;
  }

  /*
   * ========================================================================================
   * ----------------------- PSR CPI header info --------------------------------------------
   * ========================================================================================
   */
  /**
   * Get mean noise value of an item
   *
   * @param fname  Name of psr file
   * @param item   Number of item (starting with 0)
   * @param value  Mean spectral noise sample [adu]
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getItemNoise(char *fname, int item, float *value)
  {
    if (psr_openFile(fname) < 0)
      return -1;
    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d)\r\n",
              __FILE__, __func__, item);
      return -1;
    }
    if (psrfile.header.reduced) {
      *value = psrfile.item_table[item].noise;
    } else {
      fprintf(stderr, "%s:%s(): WARNING: No noise values in file\r\n",
              __FILE__, __func__);
      return -1;
    }
    return 0;
  }

  /*
   * ========================================================================================
   * ----------------------- PSR data -------------------------------------------------------
   * ========================================================================================
   */

  /**
   * Get power spectrum
   *
   * @param fname  Name of psr file
   * @param item   Item number (starting with 0)
   * @param rangegate Range gate selction (starting with 0)
   * @param complex Pointer to complex power spectrum to save
   *
   * @return 0 on success, -1 otherwise
   */
  long psr_getPowerSpectrum(char *fname, int item, int rangegate, float *complex)
  {
    void *pcpi;
    void *psrdata;
    long nrgates = -1;
    long ntxpulses = -1;

    if (psr_openFile(fname) < 0)
      return -1;

    if (item >= psrfile.header.itemtypepsr) {
      fprintf(stderr, "%s:%s(): WARNING: Itemnumber too large (%d must be smaller than %ld)\r\n",
              __FILE__, __func__, item, psrfile.header.itemtypepsr);
      return -1;
    }

    pcpi = getCpiBase(&psrfile, item);

    getValueFromCpiHeader(&psrfile, item, CPI_NUM_RANGEGATES_OFFSET, &nrgates);
    if (rangegate >= nrgates) {
      fprintf(stderr, "%s:%s(): WARNING: Rangenumber too large (%d must be smaller than %ld)\r\n",
              __FILE__, __func__, rangegate, nrgates);
      return -1;
    }
    getValueFromCpiHeader(&psrfile, item, CPI_NUM_TXPULSES_OFFSET, &ntxpulses);

    psrdata = pcpi + CPI_HEADER_LEN + rangegate*ntxpulses*COMPLEXSIZE;

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    memcpy(complex, psrdata, ntxpulses*COMPLEXSIZE);
#elif __BYTE_ORDER__ ==  __ORDER_BIG_ENDIAN__
    {
      int k;
      void *ptr = (void *)complex;
      for (k=0; k<ntxpulses*2; k++) {
        *((uint32_t *)(ptr+k*FLOATSIZE)) = le32toh(*((uint32_t *)(psrdata+k*FLOATSIZE))) & 0xffffffff;
      }
    }
#else
#error "Endianness undefined!"
#endif

    return 0;
  }

  /*
   * ========================================================================================
   * ----------------------- Other global functions -----------------------------------------
   * ========================================================================================
   */

  /**
   * Init psr module
   * Should be called at startup.
   */
  void psr_init(void)
  {
    psr_file_clean(&psrfile);
  }

  /**
   * Clean up
   * Should be called at the end
   */
  void psr_cleanup(void)
  {
    if (psrfile.filebase) {
      munmap(psrfile.filebase, psrfile.mmap_size);
      psrfile.filebase = NULL;
    }
    if (psrfile.item_table) {
      free(psrfile.item_table);
      psrfile.item_table = NULL;
    }
    if (psrfile.compressed && psrfile.pitem) {
      free(psrfile.pitem);
      psrfile.pitem = NULL;
      psrfile.itemnum = -1;
    }
    psr_file_clean(&psrfile);
  }

  /**
   * Return pointer to header structure
   * of psr file.
   *
   * @return pointer to header or NULL if no file opened.
   */
  struct psr_file_header *psr_getHeaderStruct(void)
  {
    return &psrfile.header;
  }

  /**
   * Return pointer to file structure
   * of psr file.
   *
   * @return pointer to file struct
   */
  struct psr_file *psr_getFileStruct(void)
  {
    return &psrfile;
  }

  /**
   * Open psr file, read file header and map the memory
   *
   * @param fname  Psr filename
   *
   * @return 0 on success, -1 otherwise
   */
  int psr_openFile(char *fname)
  {
    FILE *fh;
    struct stat file_stat;
    long mmap_size;

    /* Read file header if not yet read */
    if ((psrfile.header_read > 0) && (strcmp(fname, psrfile.filename) == 0) ) {
      /* header already read: nothing to do */
      return 0;
    }
    psr_cleanup();

    /* Open psr file */
    fh = fopen(fname, "r");
    if (fh == NULL) {
      fprintf(stderr, "%s:%s(): ERROR: Cannot open file \'%s\'\r\n",
              __FILE__, __func__, fname);
      return -1;
    }

    /* Read file header */
    if (psr_header_read(fh, &psrfile.header) < 0) {
      fprintf(stderr, "%s:%s(): ERROR: Cannot read header of file \'%s\'\r\n",
              __FILE__, __func__, fname);
      fclose(fh);
      return -1;
    }

    /* Get filesize of psr file*/
    if (fstat(fileno(fh), &file_stat) < 0) {
      fprintf(stderr, "%s:%s(): ERROR reading stat from file '%s': %s\r\n",
              __FILE__, __func__, fname, strerror(errno));
      fclose(fh);
      return -1;
    }
    mmap_size = file_stat.st_size;

    if (psrfile.header.reduced) {
      /* Reduced psr file */
      psrfile.item_type = ITEM_TYPE_ID_RED;

      /* Check filesize to see if the file is compressed */
      if (file_stat.st_size < psrfile.header.datalength + psrfile.header.bodyoffset) {
        psrfile.compressed = 1;

        /* Read item table only if reduced and compressed */
        if (readItemTable(&psrfile, fh) < 0) {
          fprintf(stderr, "%s:%s(): ERROR reading item table from file '%s'\r\n",
                  __FILE__, __func__, fname);
          fclose(fh);
          return -1;
        }
      }
      else if (file_stat.st_size == psrfile.header.datalength + psrfile.header.bodyoffset) {
        psrfile.compressed = 0;
      }
      else {
        fprintf(stderr, "%s:%s(): ERROR filesize too large, file '%s'\r\n",
                __FILE__, __func__, fname);
        fclose(fh);
        return -1;
      }
    }
    else {
      if (mmap_size != psrfile.header.datalength + psrfile.header.bodyoffset)  {
        fprintf(stderr, "%s:%s(): ERROR filesize missmatch, file '%s'\r\n",
                __FILE__, __func__, fname);
        fclose(fh);
        return -1;
      }
      psrfile.item_type = ITEM_TYPE_ID;
    }

    /* Map data section to memory */
    psrfile.filebase = mmap(NULL, mmap_size,
                            PROT_READ, MAP_PRIVATE,
                            fileno(fh), 0);
    if (psrfile.filebase == (void*)-1) {
      fprintf(stderr, "%s:%s(): ERROR mmap : %s\r\n",
              __FILE__, __func__, strerror(errno));
      fclose(fh);
      return -1;
    }
    fclose(fh);
    psrfile.base = psrfile.filebase + psrfile.header.bodyoffset;

    psrfile.mmap_size = mmap_size;
    strncpy(psrfile.filename, fname, MAXFILENAME);
    psrfile.header_read = 1;
    return 0;
  }

  /**
   * Read values from all cpi headers for uncompressed files
   *
   * @param ptr        File structure pointer
   * @param azStart    Pointer where the azmith start data has to be stored.
   * @param azStop     Pointer where the azimuth stop data has to be stored.
   * @param elStart    Pointer where the elevation start data has to be stored.
   * @param elStop     Pointer where the elevation stop data has to be stored.
   * @param ntxpulses  Pointer where the number of tyxpulses data has to be stored.
   * @param prf        Pointer where the PRF data has to be stored.
   * @param nrbins     Pointer where the number of range bins data has to be stored.
   * @param txpower    Pointer where the tx power data has to be stored.
   *
   * @return 0 on success, -1 other wise
   */
  static long getValueArrays_uncompressed(struct psr_file *ptr, void *azStart, void *azStop,
                                          void *elStart, void *elStop, void *ntxpulses, void *prf,
                                          void *nrbins, void *txpower)
  {
    void *pitem;
    struct item_header item;
    int k;

    pitem = ptr->base;
    k=0;
    while(k < ptr->header.itemtypepsr) {
      /* Get item info */
      item.type   = be32toh(*((uint32_t *)pitem)) & 0xffffffff;
      item.length = be64toh(*((uint64_t *)(pitem+ITEM_LENGTH_OFFSET)));

      pitem += ITEM_HEADER_OFFSET;

      if (item.type != ptr->item_type) {
        fprintf(stderr, "%s:%s() ERROR: Unexpected item type: 0x%lx\r\n",
                __FILE__, __func__, item.type);
        pitem += item.length;
        continue;
      }
      *(uint32_t *)azStart = le32toh(*((uint32_t *)(pitem + CPI_AZIMUTHSTART_OFFSET))) & 0xffffffff;
      azStart += sizeof(uint32_t);
      *(uint32_t *)azStop = le32toh(*((uint32_t *)(pitem + CPI_AZIMUTHSTOP_OFFSET))) & 0xffffffff;
      azStop += sizeof(uint32_t);
      *(uint32_t *)elStart = le32toh(*((uint32_t *)(pitem + CPI_ELEVATIONSTART_OFFSET))) & 0xffffffff;
      elStart += sizeof(uint32_t);
      *(uint32_t *)elStop = le32toh(*((uint32_t *)(pitem + CPI_ELEVATIONSTOP_OFFSET))) & 0xffffffff;
      elStop += sizeof(uint32_t);
      *(uint32_t *)ntxpulses = le32toh(*((uint32_t *)(pitem + CPI_NUM_TXPULSES_OFFSET))) & 0xffffffff;
      ntxpulses += sizeof(uint32_t);
      *(uint32_t *)prf = le32toh(*((uint32_t *)(pitem + CPI_PRF_OFFSET))) & 0xffffffff;
      prf += sizeof(uint32_t);
      *(uint32_t *)nrbins = le32toh(*((uint32_t *)(pitem + CPI_NUM_RANGEGATES_OFFSET))) & 0xffffffff;
      nrbins += sizeof(uint32_t);
      *(uint32_t *)txpower = le32toh(*((uint32_t *)(pitem + CPI_AVGTXPOWER_OFFSET))) & 0xffffffff;
      txpower += sizeof(uint32_t);

      k++;

      pitem += item.length;
    }
    return 0;
  }

  /**
   * Read values from all cpi headers for compressed files
   *
   * @param ptr        File structure pointer
   * @param azStart    Pointer where the azmith start data has to be stored.
   * @param azStop     Pointer where the azimuth stop data has to be stored.
   * @param elStart    Pointer where the elevation start data has to be stored.
   * @param elStop     Pointer where the elevation stop data has to be stored.
   * @param ntxpulses  Pointer where the number of tyxpulses data has to be stored.
   * @param prf        Pointer where the PRF data has to be stored.
   * @param nrbins     Pointer where the number of range bins data has to be stored.
   * @param txpower    Pointer where the tx power data has to be stored.
   *
   * @return 0 on success, -1 other wise
   */
  static long getValueArrays_compressed(struct psr_file *ptr, void *azStart, void *azStop,
                                        void *elStart, void *elStop, void *ntxpulses, void *prf,
                                        void *nrbins, void *txpower)
  {
    void *pitem;
    int k;

    for (k=0; k<ptr->header.itemtypepsr; k++) {
      pitem = getCpiBase(ptr, k);
      if (pitem == NULL) {
        return -1;
      }

      *(uint32_t *)azStart = le32toh(*((uint32_t *)(pitem + CPI_AZIMUTHSTART_OFFSET))) & 0xffffffff;
      azStart += sizeof(uint32_t);
      *(uint32_t *)azStop = le32toh(*((uint32_t *)(pitem + CPI_AZIMUTHSTOP_OFFSET))) & 0xffffffff;
      azStop += sizeof(uint32_t);
      *(uint32_t *)elStart = le32toh(*((uint32_t *)(pitem + CPI_ELEVATIONSTART_OFFSET))) & 0xffffffff;
      elStart += sizeof(uint32_t);
      *(uint32_t *)elStop = le32toh(*((uint32_t *)(pitem + CPI_ELEVATIONSTOP_OFFSET))) & 0xffffffff;
      elStop += sizeof(uint32_t);
      *(uint32_t *)ntxpulses = le32toh(*((uint32_t *)(pitem + CPI_NUM_TXPULSES_OFFSET))) & 0xffffffff;
      ntxpulses += sizeof(uint32_t);
      *(uint32_t *)prf = le32toh(*((uint32_t *)(pitem + CPI_PRF_OFFSET))) & 0xffffffff;
      prf += sizeof(uint32_t);
      *(uint32_t *)nrbins = le32toh(*((uint32_t *)(pitem + CPI_NUM_RANGEGATES_OFFSET))) & 0xffffffff;
      nrbins += sizeof(uint32_t);
      *(uint32_t *)txpower = le32toh(*((uint32_t *)(pitem + CPI_AVGTXPOWER_OFFSET))) & 0xffffffff;
      txpower += sizeof(uint32_t);

    }
    return 0;
  }

  /**
   * Read values from all cpi headers
   *
   * @param fname      Name of psr file
   * @param azStart    Pointer where the azmith start data has to be stored.
   * @param azStop     Pointer where the azimuth stop data has to be stored.
   * @param elStart    Pointer where the elevation start data has to be stored.
   * @param elStop     Pointer where the elevation stop data has to be stored.
   * @param ntxpulses  Pointer where the number of tyxpulses data has to be stored.
   * @param prf        Pointer where the PRF data has to be stored.
   * @param nrbins     Pointer where the number of range bins data has to be stored.
   * @param txpower    Pointer where the tx power data has to be stored.
   *
   * @return 0 on succes, -1 other wise
   */
  long psr_getValueArrays(char *fname, float *azStart, float *azStop,
                          float *elStart, float *elStop, long *ntxpulses, float *prf,
                          long *nrbins, float *txpower)
  {
    if (psr_openFile(fname) < 0)
      return -1;

    if (psrfile.compressed)
      return getValueArrays_compressed(&psrfile, azStart, azStop, elStart, elStop, ntxpulses,
                                       prf, nrbins, txpower);
    else
      return getValueArrays_uncompressed(&psrfile, azStart, azStop, elStart, elStop, ntxpulses,
                                         prf, nrbins, txpower);
  }

  /*
   * ========================= decompress functions ============================
   */

  /* states */
#define STATE_IDLE        0
#define STATE_GETNUMZEROS 1

  static uint32_t psr_decompress_numzerowrite;
  static int psr_decompress_state;
  static int psr_decompress_shift;

  /**
   * Init static variables for decompressing.
   */
  void psr_initDecompress(void)
  {
    psr_decompress_numzerowrite = 0;
    psr_decompress_state = STATE_IDLE;
    psr_decompress_shift = 0;
  }

  /**
   * Check if decompressing is done.
   */
  int psr_checkFinished()
  {
    if (psr_decompress_state == STATE_GETNUMZEROS) {
      fprintf(stderr, "Error State missmatch: Still in STATE_GETNUMZEROS at the end of file!\r\n");
      return -1;
    }
    else
      return 0;
  }

  /**
   * Decompress buffer and write the decompressed data.
   *
   * @param buf  Pointer to compressed data
   * @param int  Length of compressed data
   * @param writefunc Functionpointer to write function
   *
   * @return 0 on success, -1 otherwise
   */
  int psr_decompress(unsigned char *buf, int cnt, int(*writefunc)(void *buf, int len) )
  {
    int k, l;
    unsigned char zero = 0x00;

    for (k=0;k<cnt;k++) {

      if (buf[k] == 0) {
        /* New zero block: read the number of zeros from the next bytes */
        if (psr_decompress_state == STATE_GETNUMZEROS) {
          fprintf(stderr, "%s:%s(): ERROR State missmatch: Already in STATE_GETNUMZEROS!\r\n",
                  __FILE__, __func__);
          return -1;
        }
        psr_decompress_state = STATE_GETNUMZEROS;
        psr_decompress_numzerowrite = 0;
        psr_decompress_shift = 0;
        continue;
      }

      if (psr_decompress_state == STATE_GETNUMZEROS) {
        /* Get the number of zeros */
        psr_decompress_numzerowrite |= (buf[k] & 0x7f) << (psr_decompress_shift*7);
        if (buf[k] & 0x80) {
          psr_decompress_shift++;
          continue;
        }

        /* Write zeros */
        for (l=0; l<psr_decompress_numzerowrite; l++) {
          if (writefunc(&zero, 1) < 0) {
            fprintf(stderr, "%s:%s(): ERROR write\r\n",
                    __FILE__, __func__);
            return -1;
          }
        }
        psr_decompress_state = STATE_IDLE;
        psr_decompress_shift = 0;
        continue;
      }

      if (writefunc(&buf[k], 1) < 0) {
        fprintf(stderr, "%s:%s(): ERROR write\r\n",
                __FILE__, __func__);
        return -1;
      }
    }
    return 0;
  }

  /*
   * ========================= static functions ============================
   */

  /**
   * Read the item table from a reduced psr file.
   *
   * @param ptr   File structure pointer
   * @param fh    File handler of reduced psr file
   *
   * @return 0 on success, -1 other wise
   */
  static int readItemTable(struct psr_file *ptr, FILE *fh)
{
  long nitems;
  char line[MAXLINESIZE];
  int k;
  int num;
  int ret;
  unsigned long offset;
  unsigned long len;
  float noise;

  nitems = ptr->header.itemtypepsr;
  ptr->item_table = malloc(nitems * sizeof(struct psr_item_table_entry));
  if (ptr->item_table == NULL) {
    fprintf(stderr, "%s:%s(): ERROR: malloc failed\r\n",
            __FILE__, __func__);
    return -1;
  }

  if (fseek(fh, ptr->header.headerlen, SEEK_SET) < 0) {
    fprintf(stderr, "%s:%s(): ERROR: fseek: %s\r\n",
            __FILE__, __func__, strerror(errno));
    free(ptr->item_table);
    ptr->item_table = NULL;
    return -1;
  }

  if (fgets(line, MAXLINESIZE, fh) == NULL) { /* read dummy line */
    fprintf(stderr, "%s:%s(): ERROR reading first line of item table\r\n",
            __FILE__, __func__);
    free(ptr->item_table);
    ptr->item_table = NULL;
    return -1;
  }

  for (k = 0; k < nitems; k++) {
    if (fgets(line, MAXLINESIZE, fh) == NULL) {
      fprintf(stderr, "%s:%s(): ERROR reading item table\r\n",
              __FILE__, __func__);
      free(ptr->item_table);
      ptr->item_table = NULL;
      return -1;
    }

    num = k;          /* default for compact format */
    noise = 0.0f;
    ret = 0;

    /* Old format with noise */
    ret = sscanf(line, "i%d of=%lx l=%lx n=%g", &num, &offset, &len, &noise);
    if (ret != 4) {
      noise = 0.0f;

      /* Old format without noise, alternate spelling */
      ret = sscanf(line, "i%d of=0x%lx len=0x%lx", &num, &offset, &len);
      if (ret != 3) {
        /* Old format without noise, short spelling */
        ret = sscanf(line, "i%d of=%lx l=%lx", &num, &offset, &len);
        if (ret != 3) {
          /* New compact format: offset len noise */
          ret = sscanf(line, "%lx %lx %g", &offset, &len, &noise);
          if (ret == 3) {
            num = k;
          } else {
            noise = 0.0f;

            /* New compact format without noise: offset len */
            ret = sscanf(line, "%lx %lx", &offset, &len);
            if (ret == 2) {
              num = k;
            } else {
              if ((k == nitems - 1) && (line[0] == '\0' || line[0] == '\n')) {
                fprintf(stderr, "%s:%s(): WARNING: Item table format error\r\n",
                        __FILE__, __func__);
                continue;
              }
              fprintf(stderr, "%s:%s(): ERROR: cannot parse item table line: %s",
                      __FILE__, __func__, line);
              free(ptr->item_table);
              ptr->item_table = NULL;
              return -1;
            }
          }
        }
      }
    }

    if (num < 0 || num >= nitems) {
      fprintf(stderr, "%s:%s(): ERROR: item index out of range: %d\r\n",
              __FILE__, __func__, num);
      free(ptr->item_table);
      ptr->item_table = NULL;
      return -1;
    }

    ptr->item_table[num].offset = offset;
    ptr->item_table[num].uncompressed_len = len;
    ptr->item_table[num].noise = noise;
  }

  return 0;
}

  /**
   * Read a value from all cpi headers for uncompressed files
   *
   * @param ptr        File structure pointer
   * @param valueoffset Offset of the value (in bytes) from the
   *        start of the cpi header. Must be a multiple of 4.
   * @param data       Pointer where the dat has to be stored.
   *
   * @return 0 on success, -1 other wise
   */
  static int getValueArrayFromCpiHeader_uncompressed(struct psr_file *ptr, int valueoffset, void *data)
  {
    void *pitem;
    struct item_header item;
    int k;

    pitem = ptr->base;
    k=0;
    while(k < ptr->header.itemtypepsr) {
      /* Get item info */
      item.type   = be32toh(*((uint32_t *)pitem)) & 0xffffffff;
      item.length = be64toh(*((uint64_t *)(pitem+ITEM_LENGTH_OFFSET)));

      pitem += ITEM_HEADER_OFFSET;

      if (item.type != ptr->item_type) {
        fprintf(stderr, "%s:%s() ERROR: Unexpected item type: 0x%lx\r\n",
                __FILE__, __func__, item.type);
        pitem += item.length;
        continue;
      }
      *(uint32_t *)data = le32toh(*((uint32_t *)(pitem + valueoffset))) & 0xffffffff;
      k++;

      data += sizeof(uint32_t);
      pitem += item.length;
    }
    return 0;
  }

  /**
   * Read a value from all cpi headers for compressed files
   *
   * @param ptr        File structure pointer
   * @param valueoffset Offset of the value (in bytes) from the
   *        start of the cpi header. Must be a multiple of 4.
   * @param data       Pointer where the dat has to be stored.
   *
   * @return 0 on success, -1 other wise
   */
  static int getValueArrayFromCpiHeader_compressed(struct psr_file *ptr, int valueoffset, void *data)
  {
    void *pitem;
    int k;

    for (k=0; k<ptr->header.itemtypepsr; k++) {
      pitem = getCpiBase(ptr, k);

      *(uint32_t *)data = le32toh(*((uint32_t *)(pitem + valueoffset))) & 0xffffffff;
      data += sizeof(uint32_t);
    }
    return 0;
  }

  /**
   * Read a value from all cpi headers
   *
   * @param ptr        File structure pointer
   * @param valueoffset Offset of the value (in bytes) from the
   *        start of the cpi header. Must be a multiple of 4.
   * @param data       Pointer where the dat has to be stored.
   *
   * @return 0 on succes, -1 other wise
   */
  static int getValueArrayFromCpiHeader(struct psr_file *ptr, int valueoffset, void *data)
  {
    if (ptr->compressed)
      return getValueArrayFromCpiHeader_compressed(ptr, valueoffset, data);
    else
      return getValueArrayFromCpiHeader_uncompressed(ptr, valueoffset, data);
  }

  /**
   * Read a value from a specific cpi header
   *
   * @param ptr        File structure pointer
   * @param item       Number of item (starting with 0)
   * @param valueoffset Offset of the value (in bytes) from the
   *        start of the cpi header. Must be a multiple of 4.
   * @param data       Pointer where the dat has to be stored.
   *
   * @return 0 on succes, -1 other wise
   */
  static int getValueFromCpiHeader(struct psr_file *ptr, int itemnum, int valueoffset, void *data)
  {
    void *pcpi;

    pcpi = getCpiBase(ptr, itemnum);
    if (!pcpi)
      return -1;
    *(unsigned long *)data = le32toh(*((unsigned long *)(pcpi + valueoffset))) & 0xffffffff;

    return 0;
  }

  /**
   * Search base address of item of an uncompressed file.
   *
   * @param ptr     Pointer to loaded file
   * @param itemnum Item number (starting with 0)
   *
   * @return pointer to base address of item, null on error
   */
  static void *getCpiBase_uncompressed(struct psr_file *ptr, int itemnum)
  {
    void *pitem;
    struct item_header item;
    int k;

    if (ptr->pitem && (ptr->itemnum >= 0) && (ptr->itemnum < itemnum)) {
      pitem = ptr->pitem;
      k = ptr->itemnum;
    }
    else {
      pitem = ptr->base;
      k=0;
    }
    ptr->pitem = NULL;
    ptr->itemnum = -1;

    while(k < ptr->header.itemtypepsr) {
      /* Get item info */
      item.type   = be32toh(*((uint32_t *)pitem)) & 0xffffffff;
      item.length = be64toh(*((uint64_t *)(pitem+ITEM_LENGTH_OFFSET)));

      if (item.type != ptr->item_type) {
        fprintf(stderr, "%s:%s() ERROR: Unexpected item type: 0x%lx\r\n",
                __FILE__, __func__, item.type);
        pitem += ITEM_HEADER_OFFSET + item.length;
        continue;
      }

      if (k < itemnum) {
        pitem += ITEM_HEADER_OFFSET + item.length;
        k++;
        continue;
      }
      else {
        /* Item found */
        ptr->pitem = pitem;
        ptr->itemnum = itemnum;
        ptr->itemsize = item.length;
        break;
      }
    }
    if (ptr->pitem)
      return (ptr->pitem + ITEM_HEADER_OFFSET);
    else
      return NULL;
  }

  /**
   * Write function
   */
  static int copyToMemory(void *buf, int len)
  {
    if (len > decompress_len - decompress_position) {
      fprintf(stderr, "%s:%s() ERROR: Decompress buffer already full\r\n",
              __FILE__, __func__);
      return -1;
    }
    memcpy(decompress_buffer+decompress_position, buf, len);
    decompress_position += len;
    return 0;
  }

  /**
   * Search base address of item of an compressed file.
   *
   * @param ptr     Pointer to loaded file
   * @param itemnum Item number (starting with 0)
   *
   * @return pointer to base address of item, null on error
   */
  static void *getCpiBase_compressed(struct psr_file *ptr, int itemnum)
  {
    long offset;
    long uncompressed_len;
    long compressed_len;
    void *compressed_ptr;
    struct item_header item;

    if (ptr->pitem && (ptr->itemnum >= 0)) {
      /* Free last pitem*/
      free(ptr->pitem);
      ptr->pitem = NULL;
      ptr->itemnum = -1;
    }

    offset = ptr->item_table[itemnum].offset;
    uncompressed_len = ptr->item_table[itemnum].uncompressed_len + ITEM_HEADER_OFFSET;

    /* Get compressed size */
    if (itemnum < ptr->header.itemtypepsr-1) {
      compressed_len = ptr->item_table[itemnum+1].offset - offset;
    }
    else {
      compressed_len = ptr->mmap_size - offset;
    }

    /* Set pointer to beginning of item */
    compressed_ptr = ptr->filebase + offset;

    if (offset+compressed_len > ptr->mmap_size) {
      fprintf(stderr, "%s:%s() ERROR: Corrupt file\r\n",
              __FILE__, __func__);
      return NULL;
    }

    /* Get item info */
    item.type   = be32toh(*((uint32_t *)compressed_ptr)) & 0xffffffff;
    //item.length = be64toh(*((uint64_t *)(pitem+ITEM_LENGTH_OFFSET)));
    if (item.type != ptr->item_type) {
      fprintf(stderr, "%s:%s() ERROR: Unexpected item type: 0x%lx\r\n",
              __FILE__, __func__, item.type);
      return NULL;
    }

    /* Allocated memory */
    ptr->pitem = malloc(uncompressed_len);
    ptr->itemnum = itemnum;

    /* Decompress buffer */
    decompress_buffer = ptr->pitem;
    decompress_len = uncompressed_len;
    decompress_position = 0;

    psr_initDecompress();
    if (psr_decompress(compressed_ptr, compressed_len, copyToMemory) < 0) {
      fprintf(stderr, "%s:%s() ERROR decompressing\r\n",
              __FILE__, __func__);
      free(ptr->pitem);
      ptr->pitem = NULL;
      ptr->itemnum = -1;
      return NULL;
    }

    if (psr_checkFinished() < 0) {
      fprintf(stderr, "%s:%s() ERROR decompressing at the end\r\n",
              __FILE__, __func__);
      free(ptr->pitem);
      ptr->pitem = NULL;
      ptr->itemnum = -1;
      return NULL;
    }
    if (decompress_position != decompress_len) {
      fprintf(stderr, "%s:%s() ERROR buffer length missmatch\r\n",
              __FILE__, __func__);
      free(ptr->pitem);
      ptr->pitem = NULL;
      ptr->itemnum = -1;
      return NULL;
    }

    return (ptr->pitem + ITEM_HEADER_OFFSET);
  }

  /**
   * Search base address of item
   *
   * @param ptr     Pointer to loaded file
   * @param itemnum Item number (starting with 0)
   *
   * @return pointer to base address of item, null on error
   */
  static void *getCpiBase(struct psr_file *ptr, int itemnum)
  {
    if ((ptr->itemnum == itemnum) && ptr->pitem) {
      /* Already there */
      return (ptr->pitem + ITEM_HEADER_OFFSET);
    }

    if (ptr->compressed)
      return getCpiBase_compressed(ptr, itemnum);
    else
      return getCpiBase_uncompressed(ptr, itemnum);
  }

  /**
   * Init file structure
   *
   * @param ptr  Pointer to file structure
   */
  static void psr_file_clean(struct psr_file *ptr)
  {
    memset(ptr, 0, sizeof(struct psr_file));
    ptr->itemnum = -1;
  }

  /**
   * Print CPI header info
   *
   * @param pcpi Pointer to start of CPI header
   */
  static void printCpiHeader(void *pcpi)
  {
    struct cpi_header cpi;
    uint32_t *longptr;

    cpi.cpi_size       = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.num_rangegates = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.num_txpulses   = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.trackingNum2D  = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.syncCounter    = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.azimuthStart);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.azimuthStop);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.elevationStart);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.elevationStop);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.polarizationIndex = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.avgTxPower);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.txCorrMode     = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    longptr = (uint32_t *)(&cpi.prf);
    *longptr           = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.prfmode        = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.prfBatchRel    = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.qualityIndicator = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.testSigInjMask = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    cpi.dataTypeMask   = le32toh(*((uint32_t *)pcpi)) & 0xffffffff;
    pcpi += 4;
    printf("  cpi size=%d num_rangegates=%d num_txpulses=%d trackingNum2D=%d\r\n" \
           "      syncCounter=%d azimuthStart=%.3f azimuthStop=%.3f\r\n" \
           "      elevationStart=%.3f elevationStop=%.3f polarizationIndex=%d\r\n" \
           "      avgTxPower=%.3f txCorrMode=%d prf=%.3f prfmode=%d\r\n" \
           "      prfBatchRel=%d qualityIndicator=%d testSigInjMask=%d dataTypeMask=%d\r\n",
           cpi.cpi_size, cpi.num_rangegates, cpi.num_txpulses, cpi.trackingNum2D, \
           cpi.syncCounter, cpi.azimuthStart, cpi.azimuthStop, \
           cpi.elevationStart, cpi.elevationStop, cpi.polarizationIndex, \
           cpi.avgTxPower, cpi.txCorrMode, cpi.prf, cpi.prfmode, \
           cpi.prfBatchRel, cpi.qualityIndicator, cpi.testSigInjMask, cpi.dataTypeMask);
  }


#ifdef __cplusplus
} /* extern "C" */
#endif
