MODULE HElem
      IMPLICIT NONE
      TYPE HElement
#ifdef __CMPLX
         COMPLEX*16  v
#else
         REAL*8  v
#endif
      END TYPE
#ifdef __CMPLX
      integer, parameter :: HElementSize=2
#else
      integer, parameter :: HElementSize=1
#endif
      integer, parameter :: HElementSizeB=8*HElementSize  ! The same in bytes
      interface operator (+)
         module procedure HElemAdd
      end interface
      interface operator (-)
         module procedure HElemSub
         module procedure HElemNeg
      end interface
      interface operator (*)
         module procedure HElemMulInt
         module procedure HElemMulReal
         module procedure HElemMulHDElem
         module procedure HElemMul
      end interface
      interface operator (**)
         module procedure HElemPow
      end interface
      interface operator (/)
         module procedure HElemDiv
         module procedure HElemDivInt
      end interface
      interface operator (.AGT.)
         module procedure HElemAGT
      end interface 
      interface operator (.AGE.)
         module procedure HElemAGE
      end interface 
      interface assignment (=)
         module procedure HElemFromVal
         module procedure HElemFromZVal
         module procedure HElemFromHDElem
      end interface
      interface EXP
         module procedure HElemExp
      end interface
      interface DCONJG
         module procedure HElemDConjg
      end interface
      interface sum
          module procedure HElemSum
      end interface
      interface DREAL
         module procedure HElemDReal
      end interface
      interface LOG
         module procedure HElemLog
      end interface
      interface SQ
         module procedure HElemSq
      end interface
      interface ABS
         module procedure HElemABS
      end interface

      TYPE HDElement
         REAL*8  v
      END TYPE
      integer, parameter :: HDElementSize=1
      interface operator (+)
         module procedure HDElemAdd
      end interface
      interface operator (-)
         module procedure HDElemSub
         module procedure HDElemNeg
      end interface
      interface operator (*)
         module procedure HDElemMul
         module procedure HDElemMulInt
         module procedure HDElemMulReal
      end interface
      interface operator (/)
         module procedure HDElemDiv
         module procedure HDElemDivInt
      end interface
      interface operator (.AGT.)
         module procedure HDElemAGT
      end interface 
      interface operator (.AGE.)
         module procedure HDElemAGE
      end interface 
      interface assignment (=)
         module procedure HDElemFromHElem
         module procedure HDElemFromVal
      end interface
      interface DREAL
         module procedure HDElemDReal
      end interface
      interface EXP
         module procedure HDElemExp
      end interface
      interface LOG
         module procedure HDElemLog
      end interface
      interface SQ
         module procedure HDElemSq
      end interface
      CONTAINS

      SUBROUTINE HElemFromVal(h,v)
         TYPE(HElement), intent(out) :: h
         REAL*8, intent(in) :: v
         h%v=v
         RETURN
      END SUBROUTINE
      SUBROUTINE HElemFromZVal(h,z)
         TYPE(HElement), intent(out) :: h
         COMPLEX*16, intent(in) ::  z
         h%v=z
         RETURN
      END SUBROUTINE
      SUBROUTINE HElemFromHDElem(h,h2)
         TYPE(HElement), intent(out) :: h
         TYPE(HDElement), intent(in) ::  h2
         h%v=h2%v
         RETURN
      END SUBROUTINE
      TYPE(HElement) FUNCTION HElemAdd(h1,h2)
         TYPE(HElement), intent(in) :: h1,h2
         HElemAdd%v=h1%v+h2%v
         RETURN
      END FUNCTION
      pure type(HElement) function HElemDConjg(h)
         type(HElement), intent(in) ::  h
#ifdef __CMPLX
         HElemDConjg%v=dconjg(h%v)
#else
         HElemDConjg%v=h%v
#endif
         RETURN
      END FUNCTION
      pure type(helement) function HElemSum(h)
        type(helement), intent(in), dimension(:) :: h
        HElemSum%v = sum(h%v)
      end function
      TYPE(HElement) FUNCTION HElemSub(h1,h2)
         TYPE(HElement), intent(in) :: h1,h2
         HElemSub%v=h1%v-h2%v
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemNeg(h1)
         TYPE(HElement), intent(in) :: h1
         HElemNeg%v=-h1%v
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemPow(h,r)
         TYPE(HElement), intent(in) :: h
         REAL*8, intent(in) :: r
         HElemPow%v=h%v**r
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemMulInt(h1,i)
         TYPE(HElement), intent(in) :: h1
         integer, intent(in) :: i
         HElemMulInt%v=h1%v*i
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemMulReal(h1,r)
         TYPE(HElement), intent(in) :: h1
         real*8, intent(in) :: r
         HElemMulReal%v=h1%v*r
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemMul(h1,h2)
         TYPE(HElement), intent(in) :: h1,h2
         HElemMul%v=h1%v*h2%v
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemMulHDElem(h1,hd1)
         TYPE(HElement), intent(in) :: h1
         TYPE(HDElement), intent(in) :: hd1
         HElemMulHDElem%v=h1%v*hd1%v
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemDiv(h1,h2)
         TYPE(HElement), intent(in) :: h1,h2
         HElemDiv%v=h1%v/h2%v
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemDivInt(h1,i2)
         TYPE(HElement), intent(in) :: h1
         INTEGER, intent(in) :: i2
         HElemDivInt%v=h1%v/i2
         RETURN
      END FUNCTION

      TYPE(HElement) FUNCTION HElemExp(h1)
         TYPE(HElement), intent(in) :: h1
         HElemExp%v=EXP(h1%v)
         RETURN
      END FUNCTION
      TYPE(HElement) FUNCTION HElemLog(h1)
         TYPE(HElement), intent(in) :: h1
         HElemLog%v=LOG(h1%v)
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HElemABS(h1)
         TYPE(HElement), intent(in) :: h1
         HElemABS%v=ABS(h1%v)
         RETURN
      END FUNCTION
      REAL*8 FUNCTION HElemSq(h1)
         TYPE(HElement), intent(in) :: h1
#ifdef __CMPLX
         HElemSq=h1%v*dconjg(h1%v)
#else
         HElemSq=h1%v**2
#endif
         RETURN
      END FUNCTION
      REAL*8 FUNCTION HElemDReal(h1)
         TYPE(HElement) h1
#ifdef __CMPLX
         HElemDReal=DREAL(h1%v)
#else
         HElemDReal=h1%v
#endif
         RETURN
      END FUNCTION

      LOGICAL FUNCTION HElemAGT(h1,r2)
         TYPE(HElement), intent(in) :: h1
         REAL*8, intent(in) :: r2
         HElemAGT=ABS(h1%v).GT.r2
         RETURN
      END FUNCTION
      LOGICAL FUNCTION HElemAGE(h1,r2)
         TYPE(HElement), intent(in) :: h1
         REAL*8, intent(in) :: r2
         HElemAGE=ABS(h1%v).GE.r2
         RETURN
      END FUNCTION
      


      TYPE(HDElement) FUNCTION HDElemAdd(h1,h2)
         TYPE(HDElement), intent(in) :: h1,h2
         HDElemAdd%v=h1%v+h2%v
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemSub(h1,h2)
         TYPE(HDElement), intent(in) :: h1,h2
         HDElemSub%v=h1%v-h2%v
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemNeg(h1)
         TYPE(HDElement), intent(in) :: h1
         HDElemNeg%v=-h1%v
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemMul(h1,h2)
         TYPE(HDElement), intent(in) :: h1,h2
         HDElemMul%v=h1%v*h2%v
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemMulInt(h1,i)
         TYPE(HDElement), intent(in) :: h1
         integer, intent(in) :: i
         HDElemMulInt%v=h1%v*i
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemMulReal(h1,r)
         TYPE(HDElement), intent(in) :: h1
         real*8, intent(in) :: r
         HDElemMulReal%v=h1%v*r
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemDiv(h1,h2)
         TYPE(HDElement), intent(in) :: h1,h2
         HDElemDiv%v=h1%v/h2%v
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemDivInt(h1,i2)
         TYPE(HDElement), intent(in) :: h1
         INTEGER, intent(in) :: i2
         HDElemDivInt%v=h1%v/i2
         RETURN
      END FUNCTION

      TYPE(HDElement) FUNCTION HDElemExp(h1)
         TYPE(HDElement), intent(in) :: h1
         HDElemExp%v=EXP(h1%v)
         RETURN
      END FUNCTION
      TYPE(HDElement) FUNCTION HDElemLog(h1)
         TYPE(HDElement), intent(in) :: h1
         HDElemLog%v=LOG(h1%v)
         RETURN
      END FUNCTION
      REAL*8 FUNCTION HDElemSq(h1)
         TYPE(HDElement), intent(in) :: h1
         HDElemSq=h1%v**2
         RETURN
      END FUNCTION

      LOGICAL FUNCTION HDElemAGT(h1,r2)
         TYPE(HDElement), intent(in) :: h1
         REAL*8, intent(in) :: r2
         HDElemAGT=ABS(h1%v).GT.r2
         RETURN
      END FUNCTION
      LOGICAL FUNCTION HDElemAGE(h1,r2)
         TYPE(HDElement), intent(in) :: h1
         REAL*8, intent(in) :: r2
         HDElemAGE=ABS(h1%v).GE.r2
         RETURN
      END FUNCTION
      REAL*8 FUNCTION HDElemDReal(h)
         TYPE(HDElement) h
         HDElemDReal=h%v
         RETURN
      END FUNCTION


!  We always initialize from reals
       SUBROUTINE HDElemFromVal(h,v)
         TYPE(HDElement),intent(out) :: h
         REAL*8,intent(in) ::  v
         h%v=v
         RETURN
      END SUBROUTINE
!.. Or HElements
       SUBROUTINE HDElemFromHElem(h,h2)
         TYPE(HDElement),intent(out) :: h
         TYPE(HElement),intent(in) :: h2
#ifdef __CMPLX
! JSS This tolerance is causing problems for now.
!        IF(ABS(DIMAG(h2%v)).lt.1.D-10) THEN
!           WRITE(6,*) "Conversion from complex to real:",h2%v
!           STOP "Imaginary part non-zero"
!        ENDIF
         h%v=DREAL(h2%v)
#else
         h%v=h2%v
#endif
         RETURN
      END SUBROUTINE
END MODULE HElem
