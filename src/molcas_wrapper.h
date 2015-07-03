! This file contains all of the relevant compiler defines for _MOLCAS_ so that we
! don't have to pollute their build system

#ifdef _MOLCAS_


#define HElement_t real(dp)
#define __DOUBLERUN
#define DISABLE_FFTW
#define POINTER8
#define __INT64
#define HAVE_SSE2
#define __Linux
#define DSFMT_MEXP=19937


! The Lemony dude needs to check this
#ifdef _MOLCAS_MPP
#define PARALLEL
#define __SHARED_MEM
#endif


#ifdef _INTEL_
#define __IFORT
#endif


#endif
