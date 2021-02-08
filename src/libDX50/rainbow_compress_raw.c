/*
  *************************************************************

  Compress data to put in the BLOB of a rainbow raw data dfile

  *************************************************************

  Filename:         rainbow_compress_raw.c
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "qCompress.h"

#ifdef __cplusplus
extern "C" 
{
#endif
	/**
	* get expected compressed data size
	*
	* @param  nbytesudata	number of bytes of uncompressed data
	* @param  nbytescdata   variable where to store the number of bytes of compressed data
	*
	* @return 0 on success, -1 otherwise
	*/
	int rainbow_getCompressedDataSize(const unsigned long nbytesudata, unsigned long *nbytescdata)
	{
		return qCompressBound(nbytesudata, nbytescdata);
	}
	
	/**
	* Compress data
	*
	* @param  nbytesudata	number of bytes of uncompressed data
	* @param  udata			pointer to uncompressed data
	* @param  nbytescdata   variable where to store the number of bytes of compressed data
	* @param  cdata      	Pointer where to store the compressed data	
	*
	* @return 0 on success, -1 otherwise
	*/
	int rainbow_compressData(const unsigned long nbytesudata, const unsigned char *udata, unsigned long *nbytescdata, unsigned char *cdata)
	{	
		int i, ok;
		unsigned char* buf;
		unsigned long nbytesalloc;
		
		nbytesalloc=*nbytescdata;		
		
		ok=qCompress(nbytesudata, udata, nbytescdata, &buf);
		if (ok < 0)
		{
			printf("%s:%s: ERROR: Unable to compress data\r\n", __FILE__, __func__);
			return -1;
		}
		if (*nbytescdata > nbytesalloc)
		{
			printf("%s:%s: Size of compressed data %lu too large for allocated memory %lu\r\n", __FILE__, __func__, *nbytescdata, nbytesalloc);
			return -2;
		}
		for (i=0; i<*nbytescdata; i++)
		{ 
			//printf("%u\r\n", *cdata);
			*cdata=*buf;
			buf++;
			cdata++;
		}
		return ok; 
	}
	
	

#ifdef __cplusplus
} /* extern "C" */
#endif
