/*
  *************************************************************

  Qt uncompress function mapping

  *************************************************************

  Filename:         qUncompress.h 
  Author:           Andreas Leuenberger
  Creation date:    2012-10-25
  Last update:      2012-11-12
    
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


#ifndef __QUNCOMPRESS_H
#define __QUNCOMPRESS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned long ulong;
typedef unsigned char uchar;

int qUncompress(const unsigned char* data, int nbytes, unsigned char **dest, unsigned int *destBytes);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __QUNCOMPRESS_H */
