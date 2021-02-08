/*
  *************************************************************

  Power Spectrum Recording, Reduce and compress data

  *************************************************************

  Filename:         psr_reduce.c
  Author:           Andreas Leuenberger
  Creation date:    2012-11-23
  Last update:      2013-02-28

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

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "psr.h"

__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");

#ifdef __cplusplus
   extern "C" {
#endif

#define APPNAME "psrReduce"

#define PSRFILE_EXTENSION ".psr"
#define PSRCOMPFILE_EXTENSION ".rd"

#define MAXFNAMELEN 1024
#define CPIBUFSIZE 4096

#define SPACESPERITEM 40
#define VERSIONSTRLEN 50

/* Defines for threshold calculation: */
/**
 * We assume the squares of the spectrum samples to be exponential
 * distributed:
 *  (1) p_x(x) = 1/m exp(-x/m)
 * where m is the mean of the exponential distributed samples.
 *
 * Too determine the threshold such that about V_T (e.g. 0.98) of the noise
 * values are set to zero, we have to solve the following equation:
 *  (2) V_T = 1 - exp(-T/m)
 * where V_T is the percentage of the noise values below the threshold T.
 * The right half of the equation is the cumulative density function
 * solved for the threshold T: cdf(T) = \int_0^T p_x(x) dx = 1 - exp(-T/m)
 *
 * The mean value is equal to the total noise power P_N divided by the
 * number of pulses (= number of samples) in the spectrum.
 *  (3) m = P_N/N_pulses
 *
 * Solving the equations (2) and (3) for the threshold T:
 *  (4) T = m * ln(1/(1 - V_T)) = P_N/N_pulses * log(1/(1 - V_T))
 *
 */
#define PNOISE 50.0   /* Total noise power of whole spectrum (P_N) [adu]  */
#define PFACTOR 4.0   /* Corresponds to log(1/(1 - V_T)) in Eq. (4).
                       * A value of 4.0 corresponds to about 98% V_T */
#define PTHRESHOLD (PNOISE*PFACTOR) /* threshold = PTHRESHOLD/N_pulses */

/* Defines for noise estimation: histogram, ... */
#define NHISTOGRAM 250             /* Number of bins for histogram */
#define XMAX       5.0             /* Maximal value for histogram */
#define MIN_SAMPLES 30             /* Minimal number of tx pulses in CPI for noise calculation */
#define MIN_HISTBINS 20            /* Minimal number of histogram bins used for noise calculation */
#define DX         (XMAX/NHISTOGRAM) /* Histogram resolution */
#define MIN_RBIN  200              /* Start range bin for histogram */

static struct psr_stat {
  unsigned long items;
  unsigned long nvals;
  unsigned long nvals_zero;
  unsigned long zeros_one;
  unsigned long zeros_group[4];
  unsigned long psrfile_size;  /* [Bytes] */
  unsigned long reduced_size;  /* [Bytes] */
} psr_stat;

static unsigned long zerocnt = 0;

struct itemIndex {
  int num;
  unsigned long offset;
  unsigned long length;
  double noise;
};

/*
 * Write compacted data to file.
 *
 * @param fd  File descriptor
 */
static int writePsrFlush(int fd)
{
  int ret=0;
  unsigned char zbuf[5];
  int zlen;

  if (zerocnt > 0) {

    /* Statistics: */
    if (zerocnt == 1)                   /* stat */
      psr_stat.zeros_one++;             /* stat */

    /* write compact zero block: */
    zbuf[0] = 0;
    zlen = 1;
    while ((zerocnt > 0) && (zlen <= 4) ) {
      if (zlen > 1)
	zbuf[zlen-1] |= 0x80;
      zbuf[zlen] = zerocnt & 0x7f;
      zerocnt = zerocnt >> 7;
      zlen++;
    }

    psr_stat.zeros_group[zlen-2]++;     /* stat */
    ret = write(fd, zbuf, zlen);
  }
  return ret;
}

/**
 * Write values to file
 * Compact successive zeros
 */
static int writePsrCompress(int fd, unsigned char *buf, int len)
{
  int k;
  int ret=0;

  for (k=0; k<len;k++) {
    if (buf[k] == 0) {
      zerocnt++;
    }
    else {
      /* Write zeros */
      if (writePsrFlush(fd) < 0) {
	ret = -1;
	break;
      }
      /* write value */
      if (write(fd, buf+k, 1) < 0) {
	ret = -1;
	break;
      }
      zerocnt = 0;
    }
  }
  return ret;
}

/**
 * Write consecutive zeros to file
 *
 * @param fd : file numerator
 * @param num : Number of zeros to write
 *
 * @return 0 on success, -1 otherwise
 */
static int writeZeros(int fd, int num)
{
  if (num < 0) {
    fprintf(stderr, "%s ERROR: Negativ number of zeros to write\r\n", __func__);
    return -1;
  }
  zerocnt += num;
  return 0;
}

/**
 * Write consecutive zeros to file
 *
 * @param fd : file numerator
 * @param num : Number of zeros to write
 *
 * @return 0 on success, -1 otherwise
 */
static int fmemset(int fd, unsigned char c, int num)
{
  int k;
  int ret=0;

  if (num < 0) {
    fprintf(stderr, "%s ERROR: Negativ number to write\r\n", __func__);
    return -1;
  }
  for (k=0;k<num;k++) {
    if (writePsrCompress(fd, &c, 1) < 0) {
      ret = -1;
      break;
    }
  }
  return ret;
}

/**
 * Print usage of application
 */
static void usage(void)
{
  psr_printVersion();
  printf("Usage:\r\n");
  printf("%s <filename> [-o <outputfile>]\r\n\r\n", APPNAME);
  printf(" Reduces the data of a raw psr file. The complex spectrum values smaller than\r\n");
  printf(" the threshold are set to zero. The threshold depends on the number of samples\r\n");
  printf(" of a spectrum. In a second step, groups of consequtive zeros\r\n");
  printf(" are compressed. The reduced and compressed data are written to the outputfile.\r\n\r\n");
  printf("  <filename>     : Input file. A raw psr file.\r\n");
  printf("  -o <outputfile>: Filename where the data will be written to. If not specified,\r\n");
  printf("                   the extension '.rd' is appended to the <filename>.\r\n\r\n");
}

/**
 * psr reduce main function
 *
 * Reduce a psr file.
 *
 * @param input file
 * @param threshold
 * @param output file name
 *
 * @return 1 on success, -1 otherwise
 */
int main( int argc, char* argv[] )
{
  int opt;
  char *psrfile;
  FILE *fh;
  struct stat file_stat;

  float val_threshold_sqr;
  int outfile_set = 0;
  char outfile[MAXFNAMELEN];
  int fout;
  struct psr_file *psrin;
  int k;
  int fnamelen;
  long headerlen;
  long nitems;
  int nspaces;
  struct itemIndex *itemInfo;
  char itemstr[SPACESPERITEM];
  char versionstring[VERSIONSTRLEN];
  long item_type_red = htobe32(ITEM_TYPE_ID_RED);

  /* init statistic */
  memset(&psr_stat, 0, sizeof(psr_stat));
  zerocnt = 0;

  /* Get options */
  while ( (opt = getopt(argc, argv, "ht:o:")) != -1) {
    switch(opt) {
    case 'o':
      strncpy(outfile, optarg, MAXFNAMELEN);
      outfile_set = 1;
      break;
    case 'h':
    default:
      usage();
      exit(EXIT_SUCCESS);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "ERROR: %s: Not enough parameters\r\n\r\n", argv[0]);
    usage();
    exit(EXIT_FAILURE);
  }
  psrfile = (char *)argv[optind];

  /* Check file extension */
  fnamelen = strlen(psrfile);
  if (strcmp(psrfile+fnamelen-sizeof(PSRFILE_EXTENSION)+1, PSRFILE_EXTENSION) != 0) {
    fprintf(stderr, "WARNING: Unexpected extension of psr-file '%s'\r\n", psrfile);
  }

  /* Check if input file exist, read file size */
  fh = fopen(psrfile, "r");
  if (fh == NULL) {
    fprintf(stderr, "ERROR open file '%s': %s\r\n", psrfile, strerror(errno));
    exit(EXIT_FAILURE);
  }
  if (fstat(fileno(fh), &file_stat) < 0) {
    fprintf(stderr, "ERROR reading stat from file '%s': %s\r\n", psrfile, strerror(errno));
    fclose(fh);
    exit(EXIT_FAILURE);
  }
  psr_stat.psrfile_size = file_stat.st_size;
  fclose(fh);

  printf("=== Reducing \'%s\'\r\n", psrfile);

  /* Make out name, if not set */
  if (outfile_set == 0)
    snprintf(outfile, MAXFNAMELEN, "%s%s", psrfile, PSRCOMPFILE_EXTENSION);

  /* Check if out file alread exists */
  if (access(outfile, F_OK) == 0) {
    fprintf(stderr, "ERROR: Output file '%s' already exists. Abort!\r\n", outfile);
    exit(EXIT_FAILURE);
  }

  /* Make out file */
  fout = open(outfile, O_WRONLY | O_CREAT,
	      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
  if (fout < 0) {
    fprintf(stderr, "ERROR: Could not create file \'%s\'\r\n", outfile);
    exit(EXIT_FAILURE);
  }
  if (flock(fout, LOCK_EX) < 0) {
    fprintf(stderr, "ERROR: Cannot lock file \'%s\'\r\n", outfile);
    close(fout);
    exit(EXIT_FAILURE);
  }
  printf("--- Writing to \'%s\'\r\n", outfile);

  /* Open source file */
  if (psr_openFile(psrfile) < 0) {
    fprintf(stderr, "ERROR opening input file \'%s\'\r\n", psrfile);
    close(fout);
    exit(EXIT_FAILURE);
  }
  psrin = psr_getFileStruct();
  nitems = psrin->header.itemtypepsr;
  headerlen = psrin->header.headerlen;

  /* Copy file header */
  if (writePsrCompress(fout, psrin->filebase, headerlen) < 0) {
    fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
    close(fout);
    psr_cleanup();
    exit(EXIT_FAILURE);
  }
  /* reserve space to write an index table for the items */
  nspaces = nitems*SPACESPERITEM;
  itemInfo = malloc(nitems*sizeof(struct itemIndex));
  if (!itemInfo) {
    fprintf(stderr, "Malloc ERROR\r\n");
    close(fout);
    psr_cleanup();
    exit(EXIT_FAILURE);
  }
  if (fmemset(fout, ' ', nspaces) < 0) {
    fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
    close(fout);
    psr_cleanup();
    free(itemInfo);
    exit(EXIT_FAILURE);
  }
  if (psrin->header.bodyoffset - headerlen - nspaces < 0) {
    fprintf(stderr, "ERROR not enough space for item table\r\n");
    close(fout);
    psr_cleanup();
    free(itemInfo);
    exit(EXIT_FAILURE);
  }
  if (writeZeros(fout, psrin->header.bodyoffset - headerlen - nspaces) < 0) {
    fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
    close(fout);
    psr_cleanup();
    free(itemInfo);
    exit(EXIT_FAILURE);
  }
  //if (writePsrCompress(fout, psrin->filebase, psrin->header.bodyoffset) < 0) {
  //  fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
  //  close(fout);
  //  psr_cleanup();
  //  exit(EXIT_FAILURE);
  //}

  /* Copy PCI after PCI */
  {
    int l;
    void *pitem;
    struct item_header item;
    unsigned long nvals;
    float re, im;
    float psqr;
    void *writePtr;
    unsigned char zero[] = {0, 0, 0, 0, 0, 0, 0, 0};
#if __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
    uint32_t tmpul;
#endif
    /* Variable for histogram */
    int doNoise;
    unsigned long hist[NHISTOGRAM];
    int hind;
    double xvec[NHISTOGRAM];
    double sum;
    double sumsq;
    unsigned long htotal;
    double sum_old;
    unsigned long htotal_old;
    int ll;
    long npulses;
    long nrbins;
    void *ptr;
    double noise;
    long startsample;
    int noise_found;

    /* Init xvec for noise calculation */
    for (ll=0; ll<NHISTOGRAM; ll++) {
      xvec[ll] = ll*DX + DX/2;
    }

    pitem = psrin->base;
    k=0;
    while(k < nitems) {
      /* Get item info */
      item.type   = be32toh(*((uint32_t *)pitem));
      item.length = be64toh(*((uint64_t *)(pitem+ITEM_LENGTH_OFFSET)));

      if (item.type != ITEM_TYPE_ID) {
	fprintf(stderr, "ERROR: Unexpected item type: 0x%lx\r\n", item.type);
	pitem += item.length;
	continue;
      }

      /* Flush previous zeros */
      if (writePsrFlush(fout) < 0) {
	fprintf(stderr, "ERROR write : %s\r\n", strerror(errno));
	close(fout);
	psr_cleanup();
	free(itemInfo);
	exit(EXIT_FAILURE);
      }

      /* File item info struct */
      itemInfo[k].offset = lseek(fout, 0, SEEK_CUR);
      itemInfo[k].length = item.length;

      if (writePsrCompress(fout, (unsigned char *)&item_type_red, ITEM_LENGTH_OFFSET) < 0) {
      	fprintf(stderr, "ERROR writing item header %d to file : %s\r\n", k, strerror(errno));
      	close(fout);
      	psr_cleanup();
      	free(itemInfo);
      	exit(EXIT_FAILURE);
      }
      /* Write Item struct and PCI struct */
      if (writePsrCompress(fout, pitem+ITEM_LENGTH_OFFSET, ITEM_HEADER_OFFSET - ITEM_LENGTH_OFFSET + CPI_HEADER_LEN) < 0) {
	fprintf(stderr, "ERROR writing item header %d to file : %s\r\n", k, strerror(errno));
	close(fout);
	psr_cleanup();
	free(itemInfo);
	exit(EXIT_FAILURE);
      }

      pitem += ITEM_HEADER_OFFSET;
      /* Get number of tx pulses and number of range gates for noise calculation */
      ptr = &npulses;
      *(unsigned long *)ptr = le32toh(*((unsigned long *)(pitem + CPI_NUM_TXPULSES_OFFSET))) & 0xffffffff;
      ptr = &nrbins;
      *(unsigned long *)ptr = le32toh(*((unsigned long *)(pitem + CPI_NUM_RANGEGATES_OFFSET))) & 0xffffffff;
      pitem += CPI_HEADER_LEN;

      /* Write data of pci */
      if  (((item.length - CPI_HEADER_LEN ) % 8) != 0) {
	fprintf(stderr, "ERROR Unexpected item %d len (not a multiple of 4)\r\n", k);
	close(fout);
	psr_cleanup();
	free(itemInfo);
	exit(EXIT_FAILURE);
      }
      nvals = (item.length - CPI_HEADER_LEN) >> 3;

      /* Do noise calculation only if the minimum number of tx pulses is exceeded */
      if (npulses >= MIN_SAMPLES) {
        doNoise = 1;
        startsample = npulses*MIN_RBIN;
        memset(hist, 0, sizeof(hist));
      } else {
        doNoise = 0;
        startsample = 0;
      }

      /* Calculate threshold according to the number of pulses: */
      val_threshold_sqr = PTHRESHOLD / npulses;

      for (l=0;l<nvals;l++) {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
	memcpy(&re, pitem, FLOATSIZE);
	memcpy(&im, pitem+FLOATSIZE, FLOATSIZE);
#else
	tmpul = le32toh(*((uint32_t *)(pitem)));
	memcpy(&re, &tmpul, FLOATSIZE);
	tmpul = le32toh(*((uint32_t *)(pitem+FLOATSIZE)));
  	memcpy(&im, &tmpul, FLOATSIZE);
#endif

        psqr = re*re + im*im;
	if (psqr < val_threshold_sqr) {
	  writePtr = &zero;
	  psr_stat.nvals_zero++;
	}
	else {
	  writePtr = pitem;
	}

	psr_stat.nvals++;
	pitem += COMPLEXSIZE;

	if (writePsrCompress(fout, writePtr,  COMPLEXSIZE) < 0) {
	  fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
	  close(fout);
	  psr_cleanup();
	  free(itemInfo);
	  exit(EXIT_FAILURE);
	}

        /* Add value to histogram: Ignore too large samples */
        if (doNoise && (l >= startsample) && (psqr < XMAX)) {
          hind = (int)floorf(psqr/DX);
          hist[hind]++;
        }
      } /* << end loop spectrum samples */

      /* Calculate mean noise from histogram */
      if (doNoise) {
        htotal = hist[0];
        sum   = xvec[0] * hist[0];
        sumsq = xvec[0]*xvec[0] * hist[0];
        noise_found = 0;
        for (ll=1; ll<NHISTOGRAM; ll++) {
          sum_old = sum;
          htotal_old = htotal;
          sum    += xvec[ll] * hist[ll];
          sumsq  += xvec[ll]*xvec[ll] * hist[ll];
          htotal += hist[ll];
          if ((ll > MIN_HISTBINS) && (htotal*sumsq >= 2.*sum*sum)) {
            /* std > mean => higher values are not noise anymore! */
            noise = sum_old/htotal_old;
            noise_found = 1;
            break;
          }
        }
        if (!noise_found) {
          noise = sum/htotal;
        }
        itemInfo[k].noise = noise;
      } else {
        itemInfo[k].noise = 0.0;
      }

      k++;
    } /* << end loop over items */
  psr_stat.items = k;
  }

  /* Write last zeros */
  if (writePsrFlush(fout) < 0) {
    fprintf(stderr, "ERROR writing last zeros to '%s': %s\r\n", outfile, strerror(errno));
  }

  /* Write item info */
  lseek(fout, headerlen, SEEK_SET);
  snprintf(versionstring, VERSIONSTRLEN, "%s\n", PSR_REDUCED_STR_VERSION);
  if (write(fout, versionstring, strlen(versionstring)) < 0) {
    fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
  }
  for (k=0; k<psr_stat.items; k++) {
    snprintf(itemstr, SPACESPERITEM, "i%d of=%lx l=%lx n=%.6g\n",
	     k, itemInfo[k].offset, itemInfo[k].length, itemInfo[k].noise);
    if (write(fout, itemstr, strlen(itemstr)) < 0) {
      fprintf(stderr, "ERROR writing to file : %s\r\n", strerror(errno));
      break;
    }
  }

  /* Read file size */
  if (fstat(fout, &file_stat) < 0) {
    fprintf(stderr, "ERROR reading stat from file '%s': %s\r\n", outfile, strerror(errno));
    close(fout);
    psr_cleanup();
    free(itemInfo);
    exit(EXIT_FAILURE);
  }
  psr_stat.reduced_size = file_stat.st_size;

  close(fout);
  psr_cleanup();

  /* Print statistics: */
  printf("File sizes:\r\n  orig    : %10ld B\r\n  reduced : %10ld B (%2.1f%%)\r\n",
	 psr_stat.psrfile_size, psr_stat.reduced_size,
	 (float)psr_stat.reduced_size*100/(float)psr_stat.psrfile_size );
  printf("Complex values:\r\n  total       : %10ld\r\n  set to zero : %10ld (%2.1f%%)\r\n",
	 psr_stat.nvals, psr_stat.nvals_zero,
	 (float)psr_stat.nvals_zero*100/(float)psr_stat.nvals);
  // printf("Items       : %ld\r\n", psr_stat.items);
  printf("Number of successive zero bytes:\r\n");
  printf(  "   1     : %ld\r\n", psr_stat.zeros_one);
  for (k=0; k<4;k++) {
    printf("  %2d bit : %ld\r\n", (k+1)*7, psr_stat.zeros_group[k]);
  }

  free(itemInfo);
  exit(EXIT_SUCCESS);
}

#ifdef __cplusplus
   } /* extern "C" */
#endif
