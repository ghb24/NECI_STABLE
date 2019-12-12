#ifndef __NECI_TYPES_H
#define __NECI_TYPEs_H

#include <stdint.h>

// By default we use double precision (8-bytes per real)
typedef double real_t;

// Are we using 32- or 64-bit integers?
#ifdef INT64_
typedef int64_t int_t;
typedef uint64_t uint_t;
#else
typedef int32_t int_t;
typedef uint32_t uint_t;
#endif

class complex_sp_t {
public:
	float real;
	float imag;
};

class complex_dp_t {
public:
	double real;
	double imag;
};

class basisfn_t {
public:
	int32_t k[3];
	int32_t Ms;
	int32_t Ml;
	int32_t spacer;
	int64_t sym;
};

// A type(helement) equivalent
#ifdef CMPLX_
typedef complex_dp_t helement_t;
#else
typedef real_t helement_t;
#endif

// Just in case ...
extern "C" void stop_all (const char *, const char*);

#endif // __NECI_TYPES_H
