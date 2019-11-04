/**
 * @file dSFMT.h
 *
 * @brief double precision SIMD oriented Fast Mersenne Twister(dSFMT)
 * pseudorandom number generator based on IEEE 754 format.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University. All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#ifndef DSFMT_H
#define DSFMT_H

#include <stdio.h>
#include <assert.h>

#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
  #warning "DSFMT_MEXP is not defined. I assume DSFMT_MEXP is 19937."
#endif
  #define DSFMT_MEXP 19937
#endif
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^DSFMT_MEXP-1.
 * #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define DSFMT_N32 (DSFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#if !defined(DSFMT_BIG_ENDIAN)
#  if defined(__BYTE_ORDER) && defined(__BIG_ENDIAN)
#    if __BYTE_ORDER == __BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(_BYTE_ORDER) && defined(_BIG_ENDIAN)
#    if _BYTE_ORDER == _BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BYTE_ORDER__) && defined(__BIG_ENDIAN__)
#    if __BYTE_ORDER__ == __BIG_ENDIAN__
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(BYTE_ORDER) && defined(BIG_ENDIAN)
#    if BYTE_ORDER == BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN) \
    || defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)
#      define DSFMT_BIG_ENDIAN 1
#  endif
#endif

#if defined(DSFMT_BIG_ENDIAN) && defined(__amd64)
#  undef DSFMT_BIG_ENDIAN
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#  include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
#  if !defined(DSFMT_UINT32_DEFINED) && !defined(SFMT_UINT32_DEFINED)
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#    define UINT64_C(v) (v ## ui64)
#    define DSFMT_UINT32_DEFINED
#    if !defined(inline)
#      define inline __inline
#    endif
#  endif
#else
#  include <inttypes.h>
#  include <stdint.h>
#  if !defined(inline)
#    if defined(__GNUC__)
#      define inline __inline__
#    else
#      define inline
#    endif
#  endif
#endif

#ifndef PRIu64
#  if defined(_MSC_VER) || defined(__BORLANDC__)
#    define PRIu64 "I64u"
#    define PRIx64 "I64x"
#  else
#    define PRIu64 "llu"
#    define PRIx64 "llx"
#  endif
#endif

#ifndef UINT64_C
#  define UINT64_C(v) (v ## ULL)
#endif

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#  if !defined(__APPLE__)
#    include <altivec.h>
#  endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};

#elif defined(HAVE_SSE2)
#  include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#else  /* standard C */
/** 128-bit data structure */
union W128_T {
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/** the 128-bit internal state array */
struct DSFMT_T {
    w128_t status[DSFMT_N + 1];
    int idx;
};
typedef struct DSFMT_T dsfmt_t;

/** dsfmt internal state vector */
static dsfmt_t dsfmt_global_data;
/** dsfmt mexp for check */
extern const int dsfmt_global_mexp;

void dsfmt_gen_rand_all(dsfmt_t *dsfmt);
void dsfmt_fill_array_open_close(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_close_open(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_open_open(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_fill_array_close1_open2(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_chk_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed, int mexp);
void dsfmt_chk_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
			     int key_length, int mexp);
const char *dsfmt_get_idstring(void);
int dsfmt_get_min_array_size(void);

#if defined(__GNUC__)
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#  define DSFMT_PRE_INLINE __forceinline static
#  define DSFMT_PST_INLINE
#else
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE
#endif
DSFMT_PRE_INLINE uint32_t dsfmt_genrand_uint32(dsfmt_t *dsfmt) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_close(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE uint32_t dsfmt_gv_genrand_uint32(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_close1_open2(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_close_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_open_close(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_gv_genrand_open_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_open_close(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_close_open(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_open_open(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_fill_array_close1_open2(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_init_gen_rand(uint32_t seed) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_gv_init_by_array(uint32_t init_key[],
					     int key_length) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
					  int key_length) DSFMT_PST_INLINE;

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static uint32_t dsfmt_genrand_uint32(dsfmt_t *dsfmt) {
    uint32_t r;
    uint64_t *psfmt64 = &dsfmt->status[0].u[0];

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++] & 0xffffffffU;
    return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).  This is
 * the primitive and faster than generating numbers in other ranges.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt) {
    double r;
    double *psfmt64 = &dsfmt->status[0].d[0];

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++];
    return r;
}

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function.  This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static uint32_t dsfmt_gv_genrand_uint32(void) {
    return dsfmt_genrand_uint32(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_close1_open2(void) {
    return dsfmt_genrand_close1_open2(&dsfmt_global_data);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close_open(dsfmt_t *dsfmt) {
    return dsfmt_genrand_close1_open2(dsfmt) - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_close_open(void) {
    return dsfmt_gv_genrand_close1_open2() - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_close(dsfmt_t *dsfmt) {
    return 2.0 - dsfmt_genrand_close1_open2(dsfmt);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_open_close(void) {
    return 2.0 - dsfmt_gv_genrand_close1_open2();
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_open(dsfmt_t *dsfmt) {
    double *dsfmt64 = &dsfmt->status[0].d[0];
    union {
	double d;
	uint64_t u;
    } r;

    if (dsfmt->idx >= DSFMT_N64) {
	dsfmt_gen_rand_all(dsfmt);
	dsfmt->idx = 0;
    }
    r.d = dsfmt64[dsfmt->idx++];
    r.u |= 1;
    return r.d - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_gv_init_gen_rand() or dsfmt_gv_init_by_array() must be called
 * before this function. This function uses \b global variables.
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_gv_genrand_open_open(void) {
    return dsfmt_genrand_open_open(&dsfmt_global_data);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [1, 2) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_fill_array_close1_open2() except that this function uses
 * \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_close1_open2(double array[], int size) {
    dsfmt_fill_array_close1_open2(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1] to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() and \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_open_close(double array[], int size) {
    dsfmt_fill_array_open_close(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_close_open(double array[], int size) {
    dsfmt_fill_array_close_open(&dsfmt_global_data, array, size);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1) to the
 * specified array[] by one call. This function is the same as
 * dsfmt_gv_fill_array_close1_open2() except the distribution range.
 * This function uses \b global variables.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2() \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void dsfmt_gv_fill_array_open_open(double array[], int size) {
    dsfmt_fill_array_open_open(&dsfmt_global_data, array, size);
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 */
inline static void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed) {
    dsfmt_chk_init_gen_rand(dsfmt, seed, DSFMT_MEXP);
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed. This function uses \b global variables.
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_init_gen_rand()
 */
inline static void dsfmt_gv_init_gen_rand(uint32_t seed) {
    dsfmt_init_gen_rand(&dsfmt_global_data, seed);
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * @param dsfmt dsfmt state vector
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 */
inline static void dsfmt_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
				       int key_length) {
    dsfmt_chk_init_by_array(dsfmt, init_key, key_length, DSFMT_MEXP);
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds.
 * This function uses \b global variables.
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_init_by_array()
 */
inline static void dsfmt_gv_init_by_array(uint32_t init_key[], int key_length) {
    dsfmt_init_by_array(&dsfmt_global_data, init_key, key_length);
}

#if !defined(DSFMT_DO_NOT_USE_OLD_NAMES)
DSFMT_PRE_INLINE const char *get_idstring(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE int get_min_array_size(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void init_gen_rand(uint32_t seed) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void init_by_array(uint32_t init_key[], int key_length)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_close1_open2(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_close_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_open_close(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double genrand_open_open(void) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_open_close(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_close_open(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_open_open(double array[], int size)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void fill_array_close1_open2(double array[], int size)
    DSFMT_PST_INLINE;

/**
 * This function is just the same as dsfmt_get_idstring().
 * @return id string.
 * see also \sa dsfmt_get_idstring()
 */
inline static const char *get_idstring(void) {
    return dsfmt_get_idstring();
}

/**
 * This function is just the same as dsfmt_get_min_array_size().
 * @return minimum size of array used for fill_array functions.
 * see also \sa dsfmt_get_min_array_size()
 */
inline static int get_min_array_size(void) {
    return dsfmt_get_min_array_size();
}

/**
 * This function is just the same as dsfmt_gv_init_gen_rand().
 * @param seed a 32-bit integer used as the seed.
 * see also \sa dsfmt_gv_init_gen_rand(), \sa dsfmt_init_gen_rand().
 */
inline static void init_gen_rand(uint32_t seed) {
    dsfmt_gv_init_gen_rand(seed);
}

/**
 * This function is just the same as dsfmt_gv_init_by_array().
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * see also \sa dsfmt_gv_init_by_array(), \sa dsfmt_init_by_array().
 */
inline static void init_by_array(uint32_t init_key[], int key_length) {
    dsfmt_gv_init_by_array(init_key, key_length);
}

/**
 * This function is just the same as dsfmt_gv_genrand_close1_open2().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close1_open2() \sa
 * dsfmt_gv_genrand_close1_open2()
 */
inline static double genrand_close1_open2(void) {
    return dsfmt_gv_genrand_close1_open2();
}

/**
 * This function is just the same as dsfmt_gv_genrand_close_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_close_open() \sa
 * dsfmt_gv_genrand_close_open()
 */
inline static double genrand_close_open(void) {
    return dsfmt_gv_genrand_close_open();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_close().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_close() \sa
 * dsfmt_gv_genrand_open_close()
 */
inline static double genrand_open_close(void) {
    return dsfmt_gv_genrand_open_close();
}

/**
 * This function is just the same as dsfmt_gv_genrand_open_open().
 * @return double precision floating point number.
 * see also \sa dsfmt_genrand_open_open() \sa
 * dsfmt_gv_genrand_open_open()
 */
inline static double genrand_open_open(void) {
    return dsfmt_gv_genrand_open_open();
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_close().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_close(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_open_close(double array[], int size) {
    dsfmt_gv_fill_array_open_close(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_close_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_close_open(double array[], int size) {
    dsfmt_gv_fill_array_close_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_open_open().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_gv_fill_array_open_open(), \sa
 * dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_open_open(double array[], int size) {
    dsfmt_gv_fill_array_open_open(array, size);
}

/**
 * This function is juset the same as dsfmt_gv_fill_array_close1_open2().
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa dsfmt_fill_array_close1_open2(), \sa
 * dsfmt_gv_fill_array_close1_open2()
 */
inline static void fill_array_close1_open2(double array[], int size) {
    dsfmt_gv_fill_array_close1_open2(array, size);
}
#endif /* DSFMT_DO_NOT_USE_OLD_NAMES */

#endif /* DSFMT_H */
/**
 * @file dSFMT.cpp
 * @brief double precision SIMD-oriented Fast Mersenne Twister (dSFMT)
 * based on IEEE 754 format.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2007,2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

/** Included dSFMT.h explicitly in header - ghb24 */




#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "dSFMT-params.h"
/** dsfmt mexp for check */
static const int dsfmt_mexp = DSFMT_MEXP;

/*----------------
  STATIC FUNCTIONS
  ----------------*/
inline static uint32_t ini_func1(uint32_t x);
inline static uint32_t ini_func2(uint32_t x);
inline static void gen_rand_array_c1o2(dsfmt_t *dsfmt, w128_t *array,
				       int size);
inline static void gen_rand_array_c0o1(dsfmt_t *dsfmt, w128_t *array,
				       int size);
inline static void gen_rand_array_o0c1(dsfmt_t *dsfmt, w128_t *array,
				       int size);
inline static void gen_rand_array_o0o1(dsfmt_t *dsfmt, w128_t *array,
				       int size);
inline static int idxof(int i);
static void initial_mask(dsfmt_t *dsfmt);
static void period_certification(dsfmt_t *dsfmt);

#if defined(HAVE_SSE2)
#  include <emmintrin.h>
/** mask data for sse2 */
static __m128i sse2_param_mask;
/** 1 in 64bit for sse2 */
static __m128i sse2_int_one;
/** 2.0 double for sse2 */
static __m128d sse2_double_two;
/** -1.0 double for sse2 */
static __m128d sse2_double_m_one;

static void setup_const(void);
#endif

/**
 * This function simulate a 32-bit array index overlapped to 64-bit
 * array of LITTLE ENDIAN in BIG ENDIAN machine.
 */
#if defined(DSFMT_BIG_ENDIAN)
inline static int idxof(int i) {
    return i ^ 1;
}
#else
inline static int idxof(int i) {
    return i;
}
#endif

/**
 * This function represents the recursion formula.
 * @param r output
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param lung a 128-bit part of the internal state array
 */
#if defined(HAVE_ALTIVEC)
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
				w128_t *lung) {
    const vector unsigned char sl1 = ALTI_SL1;
    const vector unsigned char sl1_perm = ALTI_SL1_PERM;
    const vector unsigned int sl1_msk = ALTI_SL1_MSK;
    const vector unsigned char sr1 = ALTI_SR;
    const vector unsigned char sr1_perm = ALTI_SR_PERM;
    const vector unsigned int sr1_msk = ALTI_SR_MSK;
    const vector unsigned char perm = ALTI_PERM;
    const vector unsigned int msk1 = ALTI_MSK;
    vector unsigned int w, x, y, z;

    z = a->s;
    w = lung->s;
    x = vec_perm(w, (vector unsigned int)perm, perm);
    y = vec_perm(z, sl1_perm, sl1_perm);
    y = vec_sll(y, sl1);
    y = vec_and(y, sl1_msk);
    w = vec_xor(x, b->s);
    w = vec_xor(w, y);
    x = vec_perm(w, (vector unsigned int)sr1_perm, sr1_perm);
    x = vec_srl(x, sr1);
    x = vec_and(x, sr1_msk);
    y = vec_and(w, msk1);
    z = vec_xor(z, y);
    r->s = vec_xor(z, x);
    lung->s = w;
}
#elif defined(HAVE_SSE2)
/**
 * This function setup some constant variables for SSE2.
 */
static void setup_const(void) {
    static int first = 1;
    if (!first) {
	return;
    }
    sse2_param_mask = _mm_set_epi32(DSFMT_MSK32_3, DSFMT_MSK32_4,
				    DSFMT_MSK32_1, DSFMT_MSK32_2);
    sse2_int_one = _mm_set_epi32(0, 1, 0, 1);
    sse2_double_two = _mm_set_pd(2.0, 2.0);
    sse2_double_m_one = _mm_set_pd(-1.0, -1.0);
    first = 0;
}

/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param d a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *u) {
    __m128i v, w, x, y, z;

    x = a->si;
    z = _mm_slli_epi64(x, DSFMT_SL1);
    y = _mm_shuffle_epi32(u->si, SSE2_SHUFF);
    z = _mm_xor_si128(z, b->si);
    y = _mm_xor_si128(y, z);

    v = _mm_srli_epi64(y, DSFMT_SR);
    w = _mm_and_si128(y, sse2_param_mask);
    v = _mm_xor_si128(v, x);
    v = _mm_xor_si128(v, w);
    r->si = v;
    u->si = y;
}
#else /* standard C */
/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param lung a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
				w128_t *lung) {
    uint64_t t0, t1, L0, L1;

    t0 = a->u[0];
    t1 = a->u[1];
    L0 = lung->u[0];
    L1 = lung->u[1];
    lung->u[0] = (t0 << DSFMT_SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
    lung->u[1] = (t1 << DSFMT_SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
    r->u[0] = (lung->u[0] >> DSFMT_SR) ^ (lung->u[0] & DSFMT_MSK1) ^ t0;
    r->u[1] = (lung->u[1] >> DSFMT_SR) ^ (lung->u[1] & DSFMT_MSK2) ^ t1;
}
#endif

#if defined(HAVE_SSE2)
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t *w) {
    w->sd = _mm_add_pd(w->sd, sse2_double_m_one);
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1].
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0c1(w128_t *w) {
    w->sd = _mm_sub_pd(sse2_double_two, w->sd);
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0o1(w128_t *w) {
    w->si = _mm_or_si128(w->si, sse2_int_one);
    w->sd = _mm_add_pd(w->sd, sse2_double_m_one);
}
#else /* standard C and altivec */
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t *w) {
    w->d[0] -= 1.0;
    w->d[1] -= 1.0;
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1].
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0c1(w128_t *w) {
    w->d[0] = 2.0 - w->d[0];
    w->d[1] = 2.0 - w->d[1];
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0o1(w128_t *w) {
    w->u[0] |= 1;
    w->u[1] |= 1;
    w->d[0] -= 1.0;
    w->d[1] -= 1.0;
}
#endif

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_c1o2(dsfmt_t *dsfmt, w128_t *array,
				       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
		 &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
	dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	dsfmt->status[j] = array[i];
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_c0o1(dsfmt_t *dsfmt, w128_t *array,
				       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
		 &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	convert_c0o1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
	dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	dsfmt->status[j] = array[i];
	convert_c0o1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
	convert_c0o1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_o0o1(dsfmt_t *dsfmt, w128_t *array,
				       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
		 &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	convert_o0o1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
	dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	dsfmt->status[j] = array[i];
	convert_o0o1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
	convert_o0o1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_o0c1(dsfmt_t *dsfmt, w128_t *array,
				       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
		 &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
	do_recursion(&array[i], &dsfmt->status[i],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	convert_o0c1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
	dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
	do_recursion(&array[i], &array[i - DSFMT_N],
		     &array[i + DSFMT_POS1 - DSFMT_N], &lung);
	dsfmt->status[j] = array[i];
	convert_o0c1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
	convert_o0c1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param x 32-bit integer
 * @return 32-bit integer
 */
static uint32_t ini_func1(uint32_t x) {
    return (x ^ (x >> 27)) * (uint32_t)1664525UL;
}

/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param x 32-bit integer
 * @return 32-bit integer
 */
static uint32_t ini_func2(uint32_t x) {
    return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
}

/**
 * This function initializes the internal state array to fit the IEEE
 * 754 format.
 * @param dsfmt dsfmt state vector.
 */
static void initial_mask(dsfmt_t *dsfmt) {
    int i;
    uint64_t *psfmt;

    psfmt = &dsfmt->status[0].u[0];
    for (i = 0; i < DSFMT_N * 2; i++) {
        psfmt[i] = (psfmt[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
    }
}

/**
 * This function certificate the period of 2^{SFMT_MEXP}-1.
 * @param dsfmt dsfmt state vector.
 */
static void period_certification(dsfmt_t *dsfmt) {
    uint64_t pcv[2] = {DSFMT_PCV1, DSFMT_PCV2};
    uint64_t tmp[2];
    uint64_t inner;
    int i;
#if (DSFMT_PCV2 & 1) != 1
    int j;
    uint64_t work;
#endif

    tmp[0] = (dsfmt->status[DSFMT_N].u[0] ^ DSFMT_FIX1);
    tmp[1] = (dsfmt->status[DSFMT_N].u[1] ^ DSFMT_FIX2);

    inner = tmp[0] & pcv[0];
    inner ^= tmp[1] & pcv[1];
    for (i = 32; i > 0; i >>= 1) {
        inner ^= inner >> i;
    }
    inner &= 1;
    /* check OK */
    if (inner == 1) {
	return;
    }
    /* check NG, and modification */
#if (DSFMT_PCV2 & 1) == 1
    dsfmt->status[DSFMT_N].u[1] ^= 1;
#else
    for (i = 1; i >= 0; i--) {
	work = 1;
	for (j = 0; j < 64; j++) {
	    if ((work & pcv[i]) != 0) {
		dsfmt->status[DSFMT_N].u[i] ^= work;
		return;
	    }
	    work = work << 1;
	}
    }
#endif
    return;
}

/*----------------
  PUBLIC FUNCTIONS
  ----------------*/
/**
 * This function returns the identification string.  The string shows
 * the Mersenne exponent, and all parameters of this generator.
 * @return id string.
 */
const char *dsfmt_get_idstring(void) {
    return DSFMT_IDSTR;
}

/**
 * This function returns the minimum size of array used for \b
 * fill_array functions.
 * @return minimum size of array used for fill_array functions.
 */
int dsfmt_get_min_array_size(void) {
    return DSFMT_N64;
}

/**
 * This function fills the internal state array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 */
void dsfmt_gen_rand_all(dsfmt_t *dsfmt) {
    int i;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&dsfmt->status[0], &dsfmt->status[0],
		 &dsfmt->status[DSFMT_POS1], &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
	do_recursion(&dsfmt->status[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
	do_recursion(&dsfmt->status[i], &dsfmt->status[i],
		     &dsfmt->status[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [1, 2) to the
 * specified array[] by one call. The number of pseudorandom numbers
 * is specified by the argument \b size, which must be at least (SFMT_MEXP
 * / 128) * 2 and a multiple of two.  The function
 * get_min_array_size() returns this minimum size.  The generation by
 * this function is much faster than the following fill_array_xxx functions.
 *
 * For initialization, init_gen_rand() or init_by_array() must be called
 * before the first call of this function. This function can not be
 * used after calling genrand_xxx functions, without initialization.
 *
 * @param dsfmt dsfmt state vector.
 * @param array an array where pseudorandom numbers are filled
 * by this function.  The pointer to the array must be "aligned"
 * (namely, must be a multiple of 16) in the SIMD version, since it
 * refers to the address of a 128-bit integer.  In the standard C
 * version, the pointer is arbitrary.
 *
 * @param size the number of 64-bit pseudorandom integers to be
 * generated.  size must be a multiple of 2, and greater than or equal
 * to (SFMT_MEXP / 128) * 2.
 *
 * @note \b memalign or \b posix_memalign is available to get aligned
 * memory. Mac OSX doesn't have these functions, but \b malloc of OSX
 * returns the pointer to the aligned memory block.
 */
void dsfmt_fill_array_close1_open2(dsfmt_t *dsfmt, double array[], int size) {
    assert(size % 2 == 0);
    assert(size >= DSFMT_N64);
    gen_rand_array_c1o2(dsfmt, (w128_t *)array, size / 2);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1] to the
 * specified array[] by one call. This function is the same as
 * fill_array_close1_open2() except the distribution range.
 *
 * @param dsfmt dsfmt state vector.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa fill_array_close1_open2()
 */
void dsfmt_fill_array_open_close(dsfmt_t *dsfmt, double array[], int size) {
    assert(size % 2 == 0);
    assert(size >= DSFMT_N64);
    gen_rand_array_o0c1(dsfmt, (w128_t *)array, size / 2);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [0, 1) to the
 * specified array[] by one call. This function is the same as
 * fill_array_close1_open2() except the distribution range.
 *
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param dsfmt dsfmt state vector.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa fill_array_close1_open2()
 */
void dsfmt_fill_array_close_open(dsfmt_t *dsfmt, double array[], int size) {
    assert(size % 2 == 0);
    assert(size >= DSFMT_N64);
    gen_rand_array_c0o1(dsfmt, (w128_t *)array, size / 2);
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range (0, 1) to the
 * specified array[] by one call. This function is the same as
 * fill_array_close1_open2() except the distribution range.
 *
 * @param dsfmt dsfmt state vector.
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa fill_array_close1_open2()
 */
void dsfmt_fill_array_open_open(dsfmt_t *dsfmt, double array[], int size) {
    assert(size % 2 == 0);
    assert(size >= DSFMT_N64);
    gen_rand_array_o0o1(dsfmt, (w128_t *)array, size / 2);
}

#if defined(__INTEL_COMPILER)
#  pragma warning(disable:981)
#endif
/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 * @param mexp caller's mersenne expornent
 */
void dsfmt_chk_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed, int mexp) {
    int i;
    uint32_t *psfmt;

    /* make sure caller program is compiled with the same MEXP */
    if (mexp != dsfmt_mexp) {
	fprintf(stderr, "DSFMT_MEXP doesn't match with dSFMT.c\n");
	exit(1);
    }
    psfmt = &dsfmt->status[0].u32[0];
    psfmt[idxof(0)] = seed;
    for (i = 1; i < (DSFMT_N + 1) * 4; i++) {
        psfmt[idxof(i)] = 1812433253UL
	    * (psfmt[idxof(i - 1)] ^ (psfmt[idxof(i - 1)] >> 30)) + i;
    }
    initial_mask(dsfmt);
    period_certification(dsfmt);
    dsfmt->idx = DSFMT_N64;
#if defined(HAVE_SSE2)
    setup_const();
#endif
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds
 * @param dsfmt dsfmt state vector.
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 * @param mexp caller's mersenne expornent
 */
void dsfmt_chk_init_by_array(dsfmt_t *dsfmt, uint32_t init_key[],
			     int key_length, int mexp) {
    int i, j, count;
    uint32_t r;
    uint32_t *psfmt32;
    int lag;
    int mid;
    int size = (DSFMT_N + 1) * 4;	/* pulmonary */

    /* make sure caller program is compiled with the same MEXP */
    if (mexp != dsfmt_mexp) {
	fprintf(stderr, "DSFMT_MEXP doesn't match with dSFMT.c\n");
	exit(1);
    }
    if (size >= 623) {
	lag = 11;
    } else if (size >= 68) {
	lag = 7;
    } else if (size >= 39) {
	lag = 5;
    } else {
	lag = 3;
    }
    mid = (size - lag) / 2;

    psfmt32 = &dsfmt->status[0].u32[0];
    memset(dsfmt->status, 0x8b, sizeof(dsfmt->status));
    if (key_length + 1 > size) {
	count = key_length + 1;
    } else {
	count = size;
    }
    r = ini_func1(psfmt32[idxof(0)] ^ psfmt32[idxof(mid % size)]
		  ^ psfmt32[idxof((size - 1) % size)]);
    psfmt32[idxof(mid % size)] += r;
    r += key_length;
    psfmt32[idxof((mid + lag) % size)] += r;
    psfmt32[idxof(0)] = r;
    count--;
    for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
	r = ini_func1(psfmt32[idxof(i)]
		      ^ psfmt32[idxof((i + mid) % size)]
		      ^ psfmt32[idxof((i + size - 1) % size)]);
	psfmt32[idxof((i + mid) % size)] += r;
	r += init_key[j] + i;
	psfmt32[idxof((i + mid + lag) % size)] += r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % size;
    }
    for (; j < count; j++) {
	r = ini_func1(psfmt32[idxof(i)]
		      ^ psfmt32[idxof((i + mid) % size)]
		      ^ psfmt32[idxof((i + size - 1) % size)]);
	psfmt32[idxof((i + mid) % size)] += r;
	r += i;
	psfmt32[idxof((i + mid + lag) % size)] += r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % size;
    }
    for (j = 0; j < size; j++) {
	r = ini_func2(psfmt32[idxof(i)]
		      + psfmt32[idxof((i + mid) % size)]
		      + psfmt32[idxof((i + size - 1) % size)]);
	psfmt32[idxof((i + mid) % size)] ^= r;
	r -= i;
	psfmt32[idxof((i + mid + lag) % size)] ^= r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % size;
    }
    initial_mask(dsfmt);
    period_certification(dsfmt);
    dsfmt->idx = DSFMT_N64;
#if defined(HAVE_SSE2)
    setup_const();
#endif
}
#if defined(__INTEL_COMPILER)
#  pragma warning(default:981)
#endif

// Wrap around the required dSFMT functions so that they're accessible from
// fortran. We use C++'s handy reference function to allow Fortran and C to
// communicate, despite the different approaches in passing arguments.
// We only expose functions as needed.
//
// JSS (with a hat-tip to discussions with AJWT)

// This creates instantiations which may be linked, where before the
// declarations are inline static, and so don't appear in the object files

#ifdef __cplusplus
extern "C"
{
#endif
    void init_gen_rand_fwrapper(uint32_t seed)
    {
        // Initialise random number generator.
        // See also the main dSFMT code for the ability to initialise using an
        // array of seeds.
        init_gen_rand(seed);
    }

    double genrand_close_open_fwrapper(void)
    {
        // Return a random number in the interval [0,1).
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        return genrand_close_open();
    }

    void fill_array_close_open_fwrapper(double array[], int size)
    {
        // Fill an array of length size with random numbers in the interval
        // [0,1).
        // This is much faster than than repeated calls to genrand_close_open_.
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        fill_array_close_open(array, size);
    }

#ifdef __cplusplus
}
#endif

/* this routine avoids compilers reporting warnings about unused inline functions */
void dummy_dSFMT (dsfmt_t *dsfmt, int i, uint32_t a[], double b[]) {
    dsfmt_genrand_close_open(dsfmt);
    dsfmt_genrand_open_close(dsfmt);
    dsfmt_gv_genrand_uint32();
    get_idstring();
    get_min_array_size();
    init_by_array(a, i);
    genrand_close1_open2();
    genrand_open_close();
    genrand_open_open();
    fill_array_open_close(b, i);
    fill_array_open_open(b, i);
    fill_array_close1_open2(b, i);
return;
}
