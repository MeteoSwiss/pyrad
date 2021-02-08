/*
  *************************************************************

  Endianess definitions

  *************************************************************

  Filename:         endianness.h 
  Author:           Andreas Leuenberger
  Creation date:    2012-10-30
  Last update:      2012-12-19
    
  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Defines for psr reading.

  *************************************************************
  =============================================================
 */

#ifndef __ENDIANNESS_H
#define __ENDIANNESS_H

#ifdef __cplusplus
   extern "C" {
#endif

#if (!defined(__linux__) && !defined(__sun__))
#  error "Not supported OS found"
#endif

#ifndef __BYTE_ORDER__
#  ifdef __BYTE_ORDER
#    define __BYTE_ORDER__ __BYTE_ORDER
#  else
#    error "Byte order not defined"
#  endif
#endif

#ifndef __ORDER_BIG_ENDIAN__
#  ifdef __BIG_ENDIAN
#    define __ORDER_BIG_ENDIAN__ __BIG_ENDIAN
#  else 
#    error "BIG ENDIAIN not defined"
#  endif
#endif

#ifndef __ORDER_LITTLE_ENDIAN__
#  ifdef __LITTLE_ENDIAN
#    define __ORDER_LITTLE_ENDIAN__ __LITTLE_ENDIAN
#  else 
#    error "LITTLE ENDIAIN not defined"
#  endif
#endif

/* Swap bytes in 16 bit value.  */
#ifndef __bswap_constant_16
#  define __bswap_constant_16(x) \
     ((unsigned short int) ((((x) >> 8) & 0xff) | (((x) & 0xff) << 8)))
#endif

#ifdef __linux__
#  if defined __GNUC__ && __GNUC__ >= 2
#    define __bswap_16(x) \
     (__extension__							\
      ({ register unsigned short int __v, __x = (unsigned short int) (x); \
	if (__builtin_constant_p (__x))					\
	  __v = __bswap_constant_16 (__x);				\
	else								\
	  __asm__ ("rorw $8, %w0"					\
		   : "=r" (__v)						\
		   : "0" (__x)						\
		   : "cc");						\
	__v; }))
#  else
/* This is better than nothing.  */
#   define __bswap_16(x) \
     (__extension__							\
      ({ register unsigned short int __x = (unsigned short int) (x);	\
	__bswap_constant_16 (__x); }))
#  endif

#else /* __sun__ */
/* This is better than nothing.  */
#  define __bswap_16(x) \
     (__extension__                                                           \
      ({ register unsigned short int __x = (unsigned short int) (x);          \
         __bswap_constant_16 (__x); }))
#endif


/* Swap bytes in 32 bit value.  */
#ifndef __bswap_constant_32
#  define __bswap_constant_32(x) \
     ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) |	\
      (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24))
#endif

#ifdef __linux__
#  if defined __GNUC__ && __GNUC__ >= 2
#    if __WORDSIZE == 64 || (defined __i486__ || defined __pentium__             \
			 || defined __pentiumpro__ || defined __pentium4__ \
			 || defined __k8__ || defined __athlon__	\
                         || defined __k6__ || defined __nocona__             \
                         || defined __core2__ || defined __geode__           \
                         || defined __amdfam10__)
/* To swap the bytes in a word the i486 processors and up provide the
   `bswap' opcode.  On i386 we have to use three instructions.  */
#      define __bswap_32(x) \
     (__extension__							\
      ({ register unsigned int __v, __x = (x);				\
	if (__builtin_constant_p (__x))					\
	  __v = __bswap_constant_32 (__x);				\
	else								\
	  __asm__ ("bswap %0" : "=r" (__v) : "0" (__x));		\
	__v; }))
#    else
#      define __bswap_32(x)                                                       \
     (__extension__							\
      ({ register unsigned int __v, __x = (x);                                \
	if (__builtin_constant_p (__x))					\
	  __v = __bswap_constant_32 (__x);				\
	else								\
	  __asm__ ("rorw $8, %w0;"					\
		   "rorl $16, %0;"					\
		   "rorw $8, %w0"					\
		   : "=r" (__v)						\
		   : "0" (__x)						\
		   : "cc");						\
	__v; }))
#    endif
#  else
#    define __bswap_32(x) \
     (__extension__							\
      ({ register unsigned int __x = (x); __bswap_constant_32 (__x); }))
#  endif

#else /* __sun__ */
#  define __bswap_32(x) \
     (__extension__                                                           \
      ({ register unsigned int __x = (x); __bswap_constant_32 (__x); }))
#endif

#ifndef __bswap_constant_64
/* Swap bytes in 64 bit value.  */
#  define __bswap_constant_64(x) \
     ((((x) & 0xff00000000000000ull) >> 56)				\
      | (((x) & 0x00ff000000000000ull) >> 40)				\
      | (((x) & 0x0000ff0000000000ull) >> 24)				\
      | (((x) & 0x000000ff00000000ull) >> 8)				\
      | (((x) & 0x00000000ff000000ull) << 8)				\
      | (((x) & 0x0000000000ff0000ull) << 24)				\
      | (((x) & 0x000000000000ff00ull) << 40)				\
      | (((x) & 0x00000000000000ffull) << 56))
#endif

#ifdef __linux__
#  if defined __GNUC__ && __GNUC__ >= 2

#    if __WORDSIZE == 64
#      define __bswap_64(x) \
          (__extension__                                                           \
	   ({ register unsigned long __v, __x = (x);			\
	     if (__builtin_constant_p (__x))				\
	       __v = __bswap_constant_64 (__x);				\
	     else							\
	       __asm__ ("bswap %q0" : "=r" (__v) : "0" (__x));		\
	     __v; }))
#    else
#      define __bswap_64(x) \
           (__extension__							\
	    ({ union { __extension__ unsigned long long int __ll;	\
		unsigned int __l[2]; } __w, __r;			\
	      if (__builtin_constant_p (x))				\
		__r.__ll = __bswap_constant_64 (x);			\
	      else							\
		{							\
		  __w.__ll = (x);					\
		  __r.__l[0] = __bswap_32 (__w.__l[1]);			\
		  __r.__l[1] = __bswap_32 (__w.__l[0]);			\
		}							\
	      __r.__ll; }))
#    endif
#  endif

#else /* __sun__ */
#  define __bswap_64(x)							\
     (__extension__							\
      ({ union { __extension__ unsigned long long int __ll;		\
	  unsigned int __l[2]; } __w, __r;				\
	if (__builtin_constant_p (x))					\
	  __r.__ll = __bswap_constant_64 (x);				\
	else								\
	  {								\
	    __w.__ll = (x);						\
	    __r.__l[0] = __bswap_32 (__w.__l[1]);			\
	    __r.__l[1] = __bswap_32 (__w.__l[0]);			\
	  }								\
	__r.__ll; }))
#endif


#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#  define htobe16(x) __bswap_16 (x)
#  define htole16(x) (x)
#  define be16toh(x) __bswap_16 (x)
#  define le16toh(x) (x)

#  define htobe32(x) __bswap_32 (x)
#  define htole32(x) (x)
#  define be32toh(x) __bswap_32 (x)
#  define le32toh(x) (x)

#  define htobe64(x) __bswap_64 (x)
#  define htole64(x) (x)
#  define be64toh(x) __bswap_64 (x)
#  define le64toh(x) (x)

#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__

#  define htobe16(x) (x)
#  define htole16(x) __bswap_16 (x)
#  define be16toh(x) (x)
#  define le16toh(x) __bswap_16 (x)

#  define htobe32(x) (x)
#  define htole32(x) __bswap_32 (x)
#  define be32toh(x) (x)
#  define le32toh(x) __bswap_32 (x)

#  define htobe64(x) (x)
#  define htole64(x) __bswap_64 (x)
#  define be64toh(x) (x)
#  define le64toh(x) __bswap_64 (x)

#else
#  error "Unexpected byte order"
#endif

#ifdef __cplusplus
   } /* extern "C" */
#endif

#endif /* __ENDIANNESS_H */
