/*
  *************************************************************

  Implement quncompress

  *************************************************************

  Filename:         qUncompress.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-30
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

#include "qUncompress.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
    \fn QByteArray qUncompress(const QByteArray &data)

    \relates QByteArray

    Uncompresses the \a data byte array and returns a new byte array
    with the uncompressed data.

    Returns an empty QByteArray if the input data was corrupt.

    This function will uncompress data compressed with qCompress()
    from this and any earlier Qt version, back to Qt 3.1 when this
    feature was added.

    \bold{Note:} If you want to use this function to uncompress external
    data that was compressed using zlib, you first need to prepend a four
    byte header to the byte array containing the data. The header must
    contain the expected length (in bytes) of the uncompressed data,
    expressed as an unsigned, big-endian, 32-bit integer.

    \sa qCompress()
*/

/*! \relates QByteArray

    \overload

    Uncompresses the first \a nbytes of \a data and returns a new byte
    array with the uncompressed data.
*/
int qUncompress(const unsigned char* data, int nbytes, unsigned char **dest, unsigned int *destBytes)
{
    if (!data) {
        printf("qUncompress: Data is null\r\n");
        return -1;
    }
    if (nbytes <= 4) {
		if (nbytes < 4 || (data[0]!=0 || data[1]!=0 || data[2]!=0 || data[3]!=0)) {
            printf("qUncompress: Input data is corrupted\r\n");
			return -1;
		}
		else {
			dest = NULL;
			*destBytes = 0;
			return 0;
		}
    }	
    ulong expectedSize = (data[0] << 24) | (data[1] << 16) |
                       (data[2] <<  8) | (data[3]      );
    ulong len = expectedSize;
	
	// printf("qUncompress: Expected data size: %lu\r\n", expectedSize);
	if (len==0) {
		printf("qUncompress: Error. Expected data size=0\r\n");
		return -1;
	}
		
    *dest = (uchar*)malloc(len);
    if (!*dest) {
		printf("qUncompress: could not allocate enough memory to uncompress data\r\n");
		return -1;
    }

    for (;;) {		
        ulong alloc = len;
        if (len  >= (ulong)(1 << 31) ) {
            printf("qUncompress: Input data is corrupted. Uncompressed data buffer requested dimensions too large\r\n");
            return -1;
        }

        int res = uncompress(*dest, &len,
                               (uchar*)data+4, nbytes-4);		
        switch (res) {
			case Z_OK:
				if (len != alloc) {
					if (len  >= (ulong)(1 << 31)) {
						printf("qUncompress: Z_OK: Input data is corrupted. Uncompressed data buffer requested dimensions too large\r\n");
						return -1;
					}
					*dest = (uchar*)realloc(*dest, len);
					if (!*dest) {
						printf("qUncompress: Z_OK: Could not allocate enough memory to uncompress data\r\n");
						return -1;
					}
				}
				*destBytes = len;
				return 0;

			case Z_MEM_ERROR:
				printf("qUncompress: Z_MEM_ERROR: Not enough memory\r\n");
				return -1;

			case Z_BUF_ERROR:
				len *= 2;
				continue;

			case Z_DATA_ERROR:
				printf("qUncompress: Z_DATA_ERROR: Input data is corrupted\r\n");
				return -1;
        }
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif
