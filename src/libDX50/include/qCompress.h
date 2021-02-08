/*
  *************************************************************

  Qt compress function mapping

  *************************************************************

  Filename:         qCompress.h 
  Author:           Jordi Figueras i Ventura
  Creation date:    2015-01-15
  Last update:      
    
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


#ifndef __QCOMPRESS_H
#define __QCOMPRESS_H

#ifdef __cplusplus
extern "C" 
{
#endif

	typedef unsigned long ulong;
	typedef unsigned char uchar;

	int qCompress(const unsigned long sourceLen, const unsigned char *source, unsigned long *destLen, unsigned char **dest);
	int qCompressBound(const unsigned long sourceLen, unsigned long *expectedSize);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __QUNCOMPRESS_H */
