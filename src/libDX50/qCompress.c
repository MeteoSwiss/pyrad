/*
  *************************************************************

  Implement qcompress

  *************************************************************

  Filename:         qCompress.c
  Author:           Jordi Figueras i Ventura
  Creation date:    2015-01-16
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

#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

#include "qCompress.h"

#ifdef __cplusplus
extern "C" {
#endif

int qCompressBound(const unsigned long sourceLen, unsigned long *expectedSize)
{
	unsigned long len;
	
	len=compressBound(sourceLen);
	if (len  >= (ulong)(1 << 31) ) 
	{
		printf("qCompressBound: Compressed data buffer requested dimensions too large\r\n");
        return -1;
    }
	*expectedSize=len;	
	return 0;
}

int qCompress(const unsigned long sourceLen, const unsigned char *source, unsigned long *destLen, unsigned char **dest)
{
	unsigned long len, alloc;
	unsigned long expectedSize;
	
    if (!source) 
	{
        printf("qCompress: Data is null\r\n");
        return -1;
    }
	
	expectedSize=compressBound(sourceLen);
	if (expectedSize  >= (unsigned long)(1 << 31) ) 
	{
		printf("qCompress: Compressed data buffer requested dimensions too large\r\n");
        return -1;
    }
	len=expectedSize;
    
    *dest = (unsigned char*)malloc(len);
    if (!*dest) 
	{
		printf("qCompress: could not allocate enough memory to uncompress data\r\n");
		return -1;
    }

    for (;;)
	{
		alloc = len;
		if (len  >= (unsigned long)(1 << 31) ) 
		{
			printf("qCompress: Compressed data buffer requested dimensions too large\r\n");
			return -1;
		}
		int res = compress(*dest, &len, source, sourceLen);
		
		switch (res) 
		{
			case Z_OK:
				if (len != alloc) 
				{
					if (len  >= (unsigned long)(1 << 31)) {
						printf("qCompress: Z_OK: Uncompressed data buffer requested dimensions too large\r\n");
						return -1;
					}
					*dest = (unsigned char*)realloc(*dest, len);
					if (!*dest) {
						printf("qCompress: Z_OK: Could not allocate enough memory to compress data\r\n");
						return -1;
					}
				}				
				*destLen = len;				
				return 0;
			
			case Z_MEM_ERROR:
				printf("qCompress: Z_MEM_ERROR: Not enough memory\r\n");
				return -1;
				
			case Z_BUF_ERROR:
				printf("qCompress: Z_BUF_ERROR: Output buffer too small\r\n");
				len *= 2;
				continue;				
		}		
	}
   
}

#ifdef __cplusplus
} /* extern "C" */
#endif