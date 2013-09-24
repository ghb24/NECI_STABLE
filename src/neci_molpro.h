#ifdef MOLPRO

#ifndef DSFMT_MEXP
#define DSFMT_MEXP 19937
#endif

#define POINTER8
#define __INT64
#define DISABLE_FFTW

#ifndef HElement_t
#define HElement_t real(dp)
#endif

/* Parallel compilation in molpro with MPI */
/* #if defined(_MOLCAS_MPP_) && !defined(GA_TCGMSG) && !defined(GA_TCGMSG5) - TCGMSG GA no longer supported */
#if defined(_MOLCAS_MPP_)

#define PARALLEL
#define CBINDMPI

#if !defined(__OPEN64__) && !defined(__OLD_PGI__)
/* __OLD_PGI__ refers to bugs in pgi v. 7->8: see bugzilla 3616 */
#define __SHARED_MEM
#endif

#endif

#if !defined(MOLPRO_f2003) || defined(__OLD_PGI__)
#define __ISO_C_HACK
#endif

#endif /* MOLPRO */
