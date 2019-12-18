MODULE HElem
      IMPLICIT NONE

#ifdef CMPLX_
      integer, parameter :: HElement_t_size=2
#else
      integer, parameter :: HElement_t_size=1
#endif

      integer, parameter :: HElement_t_sizeB=8*HElement_t_size  ! The same in bytes

END MODULE HElem
