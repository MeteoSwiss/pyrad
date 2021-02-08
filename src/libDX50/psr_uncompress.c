/*
  *************************************************************

  Power Spectrum Recording, decompress data

  *************************************************************

  Filename:         psr_reduce.c
  Author:           Andreas Leuenberger
  Creation date:    2012-11-23
  Last update:      2013-02-21
    
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
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "psr.h"

#ifdef __cplusplus
   extern "C" {
#endif

#define APPNAME "psrUncompress"

#define PSRCOMPFILE_EXTENSION ".rd"

#define MAXFNAMELEN 1024
#define BUFSIZE 4096

static int fout = -1;

/**
 * Print usage of application
 */
static void usage(void) 
{
  psr_printVersion();
  printf("Usage:\r\n");
  printf("%s <filename> [-o <outputfile>]\r\n\r\n", APPNAME);
  printf(" Uncompress a reduced psr file. The compressed groups of consequtive zeros are\r\n");
  printf(" decompressed and written in full length.\r\n\r\n");
  printf("  <filename>     : Input file. A reduced file with psr data.\r\n");
  printf("  -o <outputfile>: Filename where the data will be written to. If not specified,\r\n");
  printf("                   the extension '.rd' is removed from the <filename>.\r\n\r\n");
}

/**
 * Write function
 */
static int writeToOutfile(void *buf, int len)
{
  int ret = write(fout, buf, len);
  if (ret < 0) {
    fprintf(stderr, "Error write to file: %s\r\n", strerror(errno));
  }
  return ret;
}

/**
 * psr uncompress main function
 *
 * Uncompress a psr file.
 * 
 * @params
 *
 * @return 0 on success, -1 otherwise
 */
int main( int argc, char* argv[] )
{
  int opt;
  char *psrfile;
  FILE *fin;
  char outfile[MAXFNAMELEN];
  int outfile_set = 0;
  unsigned char buf[BUFSIZE];
  int cnt;
  int fnamelen;

  /* Get options */
  while ( (opt = getopt(argc, argv, "ho:")) != -1) {
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
  if (strcmp(psrfile+fnamelen-sizeof(PSRCOMPFILE_EXTENSION)+1, PSRCOMPFILE_EXTENSION) != 0) {
    fprintf(stderr, "Warning: Unexpected extension of compressed psr-file '%s'\r\n", psrfile);
    if (!outfile_set) {
      fprintf(stderr, "No output filename set! Abort!\r\n");
      exit(EXIT_FAILURE);
    }
  } 

  /* Check if file exist */
  fin = fopen(psrfile, "r");
  if (fin == NULL) {
    fprintf(stderr, "Error open file '%s': %s\r\n", psrfile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  printf("-- Uncompress psr-file \'%s\'\r\n", psrfile);

  if (!outfile_set) {
    /* Make out name */
    strncpy(outfile, psrfile, MAXFNAMELEN);
    outfile[fnamelen - sizeof(PSRCOMPFILE_EXTENSION)+1] = '\0';
  }
  printf("--- Writing to \'%s\'\r\n", outfile);

  /* Check if out file alread exists */
  if (access(outfile, F_OK) == 0) {
    fprintf(stderr, "Error: Output file '%s' already exists. Abort!\r\n", outfile);
    exit(EXIT_FAILURE);
  }

  /* Make file to write */
  fout = open(outfile, O_WRONLY | O_CREAT, 
	      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
  if (fout < 0) {
    fprintf(stderr, "ERROR: Could not create file \'%s\'\r\n", outfile);
    fclose(fin);
    exit(EXIT_FAILURE);
  }

  /* Read file and decompress */
  cnt = read(fileno(fin), buf, BUFSIZE);
  if (cnt < 0) {
    fprintf(stderr, "Error reading from file '%s': %s\r\n", psrfile, strerror(errno));
    fclose(fin);
    close(fout);
    exit(EXIT_FAILURE);
  }

  psr_initDecompress();
  while (cnt > 0) {

    /* Decompress buffer and write it to file */
    if (psr_decompress(buf, cnt, writeToOutfile) < 0) {
      fprintf(stderr, "Error decompressing from file '%s'\r\n", psrfile);
      fclose(fin);
      close(fout);
      exit(EXIT_FAILURE);
    }

    /* read next buffer */
    cnt = read(fileno(fin), buf, BUFSIZE);
    if (cnt < 0) {
      fprintf(stderr, "Error reading from file '%s': %s\r\n", psrfile, strerror(errno));
      fclose(fin);
      close(fout);
      exit(EXIT_FAILURE);
    }
  }

  fclose(fin);
  close(fout);

  if (psr_checkFinished() < 0) {
    fprintf(stderr, "Error decompressing at the end of file '%s'\r\n", psrfile);
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}

#ifdef __cplusplus
   } /* extern "C" */
#endif
