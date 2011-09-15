#ifdef MOLPRO

#ifndef DSFMT_MEXP
#define DSFMT_MEXP 19937
#endif

#define POINTER8
#define HAVE_SSE2
#define __INT64
#define DISABLE_FFTW 

#ifndef HElement_t
#define HElement_t real(dp)
#endif

#ifdef _MOLCAS_MPP_
#define __SHARED_MEM
#define PARALLEL
#endif

#endif /* MOLPRO */
