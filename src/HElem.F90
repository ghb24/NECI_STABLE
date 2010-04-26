MODULE HElem
      IMPLICIT NONE

#ifdef __CMPLX
      integer, parameter :: HElement_t_size=2
#else
      integer, parameter :: HElement_t_size=1
#endif

      integer, parameter :: HElement_t_sizeB=8*HElement_t_size  ! The same in bytes

END MODULE HElem
