/*
  *************************************************************

  Power Spectrum Recording: File header (C language code)

  *************************************************************

  Filename:         psr_header.h
  Author:           Andreas Leuenberger
  Creation date:    2012-10-25
  Last update:      2020-02-25 by Zaira Schauwecker

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Defines for the header of a psr data file.

  *************************************************************
  =============================================================
 */

#ifndef __PSR_HEADER_H
#define __PSR_HEADER_H

#include <stdio.h>

#ifdef __cplusplus
   extern "C" {
#endif

#define MAXHEADER_STRLEN 32

/* File header strings */
#define PSR_HEADER_START         "[ARCHIVE_HEADER_START]"
#define PSR_HEADER_STOP          "[ARCHIVE_HEADER_END]"
#define PSR_HEADER_BODYOFFSET    "body.offset="
#define PSR_HEADER_DATALENGTH    "data.length="
#define PSR_HEADER_DATASTARTTIME "data.starttime="
#define PSR_HEADER_DATASTOPTIME  "data.stoptime="
#define PSR_HEADER_ITEMCOUNT     "item.count="
#define PSR_HEADER_ITEMTYPEPSR   "item.type.powerspectrum="
#define PSR_HEADER_RSPPRF        "states.rspprf="
#define PSR_HEADER_RSPSTART      "states.rsprstart="
#define PSR_HEADER_RSPSTOP       "states.rsprstop="
#define PSR_HEADER_RSPSTEP       "states.rsprstep="
#define PSR_HEADER_RSPPWIDTH     "states.rsppwidth="
#define PSR_HEADER_ADU2DBMOFFSET_H "states.spbdbmtologoffset="
#define PSR_HEADER_ADU2DBMOFFSET_V "states.spbdpvdbmtologoffset="
#define PSR_HEADER_RADCONST_H    "states.rspdphradconst="
#define PSR_HEADER_RADCONST_V    "states.rspdpvradconst="
#define PSR_HEADER_NOISEPWR_H    "states.rspnoisepwr="
#define PSR_HEADER_NOISEPWR_V    "states.rspdpvnoisepwr="
#define PSR_HEADER_PATHATT       "states.rspathatt="
#define PSR_HEADER_MFLOSS        "states.gdrxmfloss="
#define PSR_HEADER_RXLOSS_V      "states.spbdpvrxloss="
#define PSR_HEADER_TXLOSS_V      "states.spbdpvtxloss="
#define PSR_HEADER_RADOMELOSS    "states.spbradomloss="
#define PSR_HEADER_RXLOSS_H      "states.spbrxloss="
#define PSR_HEADER_TXLOSS_H      "states.spbtxloss="
#define PSR_HEADER_MFBANDWIDTH   "states.gdrxmfwidth="
#define PSR_HEADER_TXPOWFAC_H    "states.spbdphrelpowtx="
#define PSR_HEADER_TXPOWFAC_V    "states.spbdpvrelpowtx="
#define PSR_HEADER_TXPOWNORM     "states.gdrxmaxpowerkwpw=" /* since SP update on October 2019 state gdrxmaxpwrkwpw is used instead of spbtxpowkw */


#define NUM_PW 4  /* number of pulse widths*/

/* File header information */
typedef struct psr_file_header {
  unsigned long bodyoffset;
  unsigned long datalength;
  char datastarttime[MAXHEADER_STRLEN];
  char datastoptime[MAXHEADER_STRLEN];
  unsigned long itemcount;
  unsigned long itemtypepsr;
  unsigned long pwidth_id;
  //unsigned long rspprf;
  float rspstart;
  float rspstop;
  float rspstep;
  float pathatt;
  float rxloss_h;
  float rxloss_v;
  float txloss_h;
  float txloss_v;
  float radomeloss;
  float txpowerfact_h;
  float txpowerfact_v;
  float adu2dbmoffset_h[NUM_PW];
  float adu2dbmoffset_v[NUM_PW];
  float radconst_h[NUM_PW];
  float radconst_v[NUM_PW];
  float noisepwr_h[NUM_PW];
  float noisepwr_v[NUM_PW];
  float mfloss[NUM_PW];
  float mfbandwidth[NUM_PW];
  float txnormpower[NUM_PW];
  long headerlen;
  unsigned char reduced;
} PSRHDR;

int psr_header_read(FILE *fh, struct psr_file_header *pheader);

#ifdef __cplusplus
   } /* extern "C" */
#endif

#endif /* __PSR_HEADER_H */
