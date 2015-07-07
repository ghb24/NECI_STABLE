#ifdef _MOLCAS_

#define HElement_t real(dp)
#define __DOUBLERUN
#define DISABLE_FFTW
#define POINTER8
#define __INT64
#define HAVE_SSE2
#define __Linux
#define DSFMT_MEXP 19937
#define DDISABLE_FFTW
#define _CONFIG Molcas

#ifdef _MOLCAS_MPP
#define PARALLEL
#endif

#ifdef _INTEL_
#define __IFORT
#endif

#endif
