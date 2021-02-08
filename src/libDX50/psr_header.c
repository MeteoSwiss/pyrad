/*
  *************************************************************

  Power Spectrum Recording: File header (C language code)

  *************************************************************

  Filename:         psr_header.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-25
  Last update:      2020-03-03 by Zaira Schauwecker (header_table line numbers adjusted)

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Read the header of a psr data file.

  *************************************************************
  =============================================================
 */

#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "psr_header.h"
#include "psr.h"

#ifdef __cplusplus
   extern "C" {
#endif

/* Header tags to read */
static const struct {
  char *searchstring;
  int   strlen;
  enum {ULONG, STRING, FLOAT, FLTARR4, END} type;
  int offset;
  unsigned int  flag;
  unsigned int linenum;
} header_table[] = {
  /* Must be ordered by linenumber! */
  /* searchstring            strlen                              type     offset                                flag    linenum */
  {PSR_HEADER_BODYOFFSET,    sizeof(PSR_HEADER_BODYOFFSET)-1,    ULONG,   offsetof(PSRHDR, bodyoffset),         1<< 0,   2},
  {PSR_HEADER_DATALENGTH,    sizeof(PSR_HEADER_DATALENGTH)-1,    ULONG,   offsetof(PSRHDR, datalength),         1<< 1,   3},
  {PSR_HEADER_DATASTARTTIME, sizeof(PSR_HEADER_DATASTARTTIME)-1, STRING,  offsetof(PSRHDR, datastarttime),      1<< 2,   4},
  {PSR_HEADER_DATASTOPTIME,  sizeof(PSR_HEADER_DATASTOPTIME)-1,  STRING,  offsetof(PSRHDR, datastoptime),       1<< 3,   5},
  {PSR_HEADER_ITEMCOUNT,     sizeof(PSR_HEADER_ITEMCOUNT)-1,     ULONG,   offsetof(PSRHDR, itemcount),          1<< 4,   6},
  {PSR_HEADER_ITEMTYPEPSR,   sizeof(PSR_HEADER_ITEMTYPEPSR)-1,   ULONG,   offsetof(PSRHDR, itemtypepsr),        1<< 5,   7},
  {PSR_HEADER_TXPOWNORM,     sizeof(PSR_HEADER_TXPOWNORM)-1,     FLTARR4, offsetof(PSRHDR, txnormpower),        1<< 6,  45}, /* since SP update on October 2019 state gdrxmaxpwrkwpw is used instead of spbtxpowkw */
  {PSR_HEADER_MFLOSS,        sizeof(PSR_HEADER_MFLOSS)-1,        FLTARR4, offsetof(PSRHDR, mfloss),             1<< 7,  54},
  {PSR_HEADER_MFBANDWIDTH,   sizeof(PSR_HEADER_MFBANDWIDTH)-1,   FLTARR4, offsetof(PSRHDR, mfbandwidth),        1<< 8,  59},
  {PSR_HEADER_PATHATT,       sizeof(PSR_HEADER_PATHATT)-1,       FLOAT,   offsetof(PSRHDR, pathatt),            1<< 9, 129},
  {PSR_HEADER_RADCONST_H,    sizeof(PSR_HEADER_RADCONST_H)-1,    FLTARR4, offsetof(PSRHDR, radconst_h),         1<<10, 135},
  {PSR_HEADER_NOISEPWR_V,    sizeof(PSR_HEADER_NOISEPWR_V)-1,    FLTARR4, offsetof(PSRHDR, noisepwr_v),         1<<11, 136},
  {PSR_HEADER_RADCONST_V,    sizeof(PSR_HEADER_RADCONST_V)-1,    FLTARR4, offsetof(PSRHDR, radconst_v),         1<<12, 139},
  {PSR_HEADER_NOISEPWR_H,    sizeof(PSR_HEADER_NOISEPWR_H)-1,    FLTARR4, offsetof(PSRHDR, noisepwr_h),         1<<13, 146},
  {PSR_HEADER_RSPPWIDTH,     sizeof(PSR_HEADER_RSPPWIDTH)-1,     ULONG,   offsetof(PSRHDR, pwidth_id),          1<<14, 160},
  {PSR_HEADER_RSPSTART,      sizeof(PSR_HEADER_RSPSTART)-1,      FLOAT,   offsetof(PSRHDR, rspstart),           1<<15, 164},
  {PSR_HEADER_RSPSTEP,       sizeof(PSR_HEADER_RSPSTEP)-1,       FLOAT,   offsetof(PSRHDR, rspstep),            1<<16, 165},
  {PSR_HEADER_RSPSTOP,       sizeof(PSR_HEADER_RSPSTOP)-1,       FLOAT,   offsetof(PSRHDR, rspstop),            1<<17, 166},
  {PSR_HEADER_ADU2DBMOFFSET_H,sizeof(PSR_HEADER_ADU2DBMOFFSET_H)-1, FLTARR4, offsetof(PSRHDR, adu2dbmoffset_h), 1<<18, 386},
  {PSR_HEADER_TXPOWFAC_H,    sizeof(PSR_HEADER_TXPOWFAC_H)-1,    FLOAT,   offsetof(PSRHDR, txpowerfact_h),      1<<19, 401},
  {PSR_HEADER_ADU2DBMOFFSET_V,sizeof(PSR_HEADER_ADU2DBMOFFSET_V)-1, FLTARR4, offsetof(PSRHDR, adu2dbmoffset_v), 1<<20, 411},
  {PSR_HEADER_TXPOWFAC_V,    sizeof(PSR_HEADER_TXPOWFAC_V)-1,    FLOAT,   offsetof(PSRHDR, txpowerfact_v),      1<<21, 431},
  {PSR_HEADER_RXLOSS_V,      sizeof(PSR_HEADER_RXLOSS_V)-1,      FLOAT,   offsetof(PSRHDR, rxloss_v),           1<<22, 435},
  {PSR_HEADER_TXLOSS_V,      sizeof(PSR_HEADER_TXLOSS_V)-1,      FLOAT,   offsetof(PSRHDR, txloss_v),           1<<23, 441},
  {PSR_HEADER_RADOMELOSS,    sizeof(PSR_HEADER_RADOMELOSS)-1,    FLOAT,   offsetof(PSRHDR, radomeloss),         1<<24, 567},
  {PSR_HEADER_RXLOSS_H,      sizeof(PSR_HEADER_RXLOSS_H)-1,      FLOAT,   offsetof(PSRHDR, rxloss_h),           1<<25, 602},
  {PSR_HEADER_TXLOSS_H,      sizeof(PSR_HEADER_TXLOSS_H)-1,      FLOAT,   offsetof(PSRHDR, txloss_h),           1<<26, 647},
  {PSR_HEADER_STOP,          sizeof(PSR_HEADER_STOP)-1,          END,     offsetof(PSRHDR, headerlen),          1<<27, 804},
};
#define NUM_HEADER_INFOS (sizeof(header_table) / sizeof(header_table[0]))
#define ALL_HEADER_INFO_FLAGS ((1<<NUM_HEADER_INFOS) -1)

/**
 * Check the file format: reduced or not
 *
 */
static unsigned char checkReduced(FILE *fh, int offset)
{
  char line[MAXLINESIZE];

  if (fseek(fh, offset, SEEK_SET) < 0) {
    fprintf(stderr, "%s:%s(): ERROR: fseek: %s\r\n",
            __FILE__, __func__, strerror(errno));
  }
  if (fgets(line, MAXLINESIZE, fh) != NULL) {
    if (strncmp(line, PSR_REDUCED_STR, sizeof(PSR_REDUCED_STR)-1) == 0) {
      return 1;
    }
  }
  return 0;
}

/**
 * Get header value
 */
static unsigned int getHeaderValue(FILE *fh, char *line, int num, void *hdr)
{
  int slen;
  float *fltarr;

  switch (header_table[num].type) {
  case ULONG:
    if (sscanf(line + header_table[num].strlen, "%d", (uint32_t *)(hdr+header_table[num].offset)) != 1) {
      fprintf(stderr, "%s:%s(): ERROR: Cannot read \'%s\'. Unkown format.\r\n",
             __FILE__, __func__, header_table[num].searchstring);
      return -1;
    }
    /* printf("--- %s %ld\r\n", header_table[num].searchstring,
     *((uint32_t *)(hdr+header_table[num].offset))); */
    break;
  case FLOAT:
    if (sscanf(line + header_table[num].strlen, "%f", (float *)(hdr+header_table[num].offset)) != 1) {
      fprintf(stderr, "%s:%s(): ERROR: Cannot read \'%s\'. Unkown format.\r\n",
             __FILE__, __func__, header_table[num].searchstring);
      return -1;
    }
    /* printf("--- %s %.3f\r\n", header_table[num].searchstring,
     *((float *)(hdr+header_table[num].offset))); */
    break;
  case FLTARR4:
    fltarr = (float *)(hdr+header_table[num].offset);
    if (sscanf(line + header_table[num].strlen, "%f %f %f %f", fltarr, fltarr+1, fltarr+2, fltarr+3) != 4) {
      fprintf(stderr, "%s:%s(): ERROR: Cannot read \'%s\'. Unkown format!\r\n",
             __FILE__, __func__, header_table[num].searchstring);
      return -1;
    }
    /* printf("--- %s %.3f\r\n", header_table[num].searchstring,
     *((float *)(hdr+header_table[num].offset))); */
    break;
  case STRING:
    /* delete \r\n at the end of the line */
    slen = strlen(line + header_table[num].strlen);
    if ((strlen > 0) && (line[header_table[num].strlen + slen - 1] == '\n'))
      line[header_table[num].strlen + slen - 1] = '\0';
    strncpy((char *)(hdr+header_table[num].offset), line + header_table[num].strlen, MAXHEADER_STRLEN);
    /* printf("--- %s %s", header_table[num].searchstring,
       (char *)(hdr+header_table[num].offset)); */
    break;
  case END: {
    /* Get the size of the header */
    long *tmplongp;
    tmplongp = (long *)(hdr+header_table[num].offset);
    *tmplongp = ftell(fh);
    break;
  }
  default:
    fprintf(stderr, "%s:%s(): ERROR: Unexpected header type %d\r\n",
           __FILE__, __func__, header_table[num].type);
    return -1;
  }
  return 1;
}

/**
 * Read psr file header
 *
 * @param fh  File handler to psr data file (position must be at the
 *            beginning of the file)
 * @param pheader Pointer to header structure
 *
 * @return 0 on success, -1 otherwise
 */
int psr_header_read(FILE *fh, struct psr_file_header *pheader)
{
  char line[MAXLINESIZE];
  int k;
  unsigned int foundflags = 0;
  unsigned int foundcnt;
  int linenumber = 0;
  int hdrstr_cnt = 0;

  /* read header */
  while(fgets(line, MAXLINESIZE, fh) != NULL) {
    linenumber++;

    if (linenumber < header_table[hdrstr_cnt].linenum)
      continue; /* skip line (no header info expected) */

    if (strncmp(line, header_table[hdrstr_cnt].searchstring, header_table[hdrstr_cnt].strlen) == 0) {
      /* string found on expected line: read value */
      getHeaderValue(fh, line, hdrstr_cnt, (void *)pheader);
    }
    else {
      fprintf(stderr, "%s:%s(): WARNING: Header string \'%s\' not found on line %d\r\n",
             __FILE__, __func__, header_table[hdrstr_cnt].searchstring, header_table[hdrstr_cnt].linenum);
      break; /* linenumber search is not working: do slow header string search. */
    }

    hdrstr_cnt++; /* next header string */
    if (hdrstr_cnt >= NUM_HEADER_INFOS) {
      /* all header infos found. Get format */
      pheader->reduced = checkReduced(fh, pheader->headerlen);
      return 0;
    }
  }

  /* ============================================================= */
  /* --------------- Do slow header read ------------------------- */
  /* ============================================================= */

  /* Not all header strings found by line number. Do it linear */
  /* Set at the beginning of the file */
  fseek(fh, 0, SEEK_SET);

  foundcnt = 0;
  linenumber = 0;
  while(fgets(line, MAXLINESIZE, fh) != NULL) {
    linenumber++;

    k=0;
    while (k<NUM_HEADER_INFOS) {
      if (strncmp(line, header_table[k].searchstring, header_table[k].strlen) == 0) {

        printf("Header string \'%s\' found on line %d\r\n", header_table[k].searchstring, linenumber);

        if (getHeaderValue(fh, line, k, (void *)pheader) > 0) {
          foundflags |= header_table[k].flag;
          foundcnt++;
          break;
        }
      }
      k++;
    }

    if (foundcnt >= NUM_HEADER_INFOS) {
      break;
    }

    if (strncmp(line, PSR_HEADER_STOP, sizeof(PSR_HEADER_STOP)-1) == 0) {
      /* printf("Header stop found\r\n"); */
      break;
    }
  }

  if (foundflags != ALL_HEADER_INFO_FLAGS) {
    fprintf(stderr, "%s:%s(): WARNING: Not all header tags found (0x%04x instead of 0x%04x)\r\n",
            __FILE__, __func__, foundflags, ALL_HEADER_INFO_FLAGS);
  }

  /* Check file format */
  pheader->reduced = checkReduced(fh, pheader->headerlen);

  return 0;
}

#ifdef __cplusplus
   } /* extern "C" */
#endif
