MODULE HElement
      IMPLICIT NONE
      TYPE HElement
#ifdef __CMPLX
         COMPLEX*16  v
#else
         REAL*8  v
#endif
      END TYPE
#ifdef __CMPLX
      PARAMETER HElementSize=2
#else
      PARAMETER HElementSize=1
#endif
      interface operator (+)
         module procedure HElemAdd
      end interface
      interface operator (-)
         module procedure HElemSub
         module procedure HElemNeg
      end interface
      interface operator (*)
         module procedure HElemMul
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
      interface LOG
         module procedure HElemLog
      end interface
      interface SQ
         module procedure HElemSq
      end interface

      TYPE HDElement
         REAL*8  v
      END TYPE
      PARAMETER HDElementSize=1
      interface operator (+)
         module procedure HDElemAdd
      end interface
      interface operator (-)
         module procedure HDElemSub
         module procedure HDElemNeg
      end interface
      interface operator (*)
         module procedure HDElemMul
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

!  We always initialize from reals
      SUBROUTINE HElemFromVal(h,v)
         TYPE(HElement) h
         REAL*8  v
         h%v=v
         RETURN
      END
      SUBROUTINE HElemFromZVal(h,z)
         TYPE(HElement) h
         COMPLEX*16  z
         h%v=z
         RETURN
      END
      SUBROUTINE HElemFromHDElem(h,h2)
         TYPE(HElement) h
         TYPE(HDElement) h2
         h%v=h2%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemAdd(h1,h2)
         TYPE(HElement) h1,h2
         HElemAdd%v=h1%v+h2%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemDConjg(h)
         TYPE(HElement) h
         HElemDConjg%v=DCONJG(h%v)
         RETURN
      END
      TYPE(HElement) FUNCTION HElemSub(h1,h2)
         TYPE(HElement) h1,h2
         HElemSub%v=h1%v-h2%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemNeg(h1)
         TYPE(HElement) h1
         HElemNeg%v=-h1%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemMul(h1,h2)
         TYPE(HElement) h1,h2
         HElemMul%v=h1%v*h2%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemDiv(h1,h2)
         TYPE(HElement) h1,h2
         HElemDiv%v=h1%v/h2%v
         RETURN
      END
      TYPE(HElement) FUNCTION HElemDivInt(h1,i2)
         TYPE(HElement) h1
         INTEGER i2
         HElemDivInt%v=h1%v/i2
         RETURN
      END

      TYPE(HElement) FUNCTION HElemExp(h1)
         TYPE(HElement) h1
         HElemExp%v=EXP(h1%v)
         RETURN
      END
      TYPE(HElement) FUNCTION HElemLog(h1)
         TYPE(HElement) h1
         HElemLog%v=LOG(h1%v)
         RETURN
      END
      REAL*8 FUNCTION HElemSq(h1)
         TYPE(HElement) h1
         HElemSq=h1%v**2
         RETURN
      END

      LOGICAL FUNCTION HElemAGT(h1,r2)
         TYPE(HElement) h1
         REAL*8 r2
         HElemAGT=ABS(h1%v).GT.r2
         RETURN
      END
      LOGICAL FUNCTION HElemAGE(h1,r2)
         TYPE(HElement) h1
         REAL*8 r2
         HElemAGE=ABS(h1%v).GE.r2
         RETURN
      END
      


!.. GETHELEMENT2
!.. Get matrix element of the hamiltonian
!.. IC is the number of basis fns that differ in NI and NJ (or -1 if not known)
!.. ECORE is the uniform background energy

      TYPE(HElement) FUNCTION GetHElement2(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,TMat,NMAX,ALAT,UMat,iC2,ECore)
         INTEGER NMSH,NMAX
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*)
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(*)
         TYPE(HElement) TMat(*),UMat(*)
         INTEGER I,nEl,NI(nEl),NJ(nEl),iC,nBasisMax(5,2),iC2
         REAL*8 ECore
         TYPE(HElement) Sum,Sum2
         INTEGER IGETEXCITLEVEL
         LOGICAL ISCSF
         INTEGER ISUB
         IF(ISCSF(NI,NEL).OR.ISCSF(NJ,NEL)) THEN
            CALL CSFGETHELEMENT(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,TMat,NMAX,ALAT,UMat,ECore,Sum2)
            GETHELEMENT2=SUM2
            RETURN
         ENDIF
         IC=IC2
         GetHElement2%v=0.D0
         IF(IC.LT.0) IC=IGETEXCITLEVEL(NI,NJ,NEL)
!.. if we differ by more than 2 spin orbital, then the hamiltonian element is 0         
         IF(IC.GT.2) RETURN
!.. SLTCND has IC is # electrons the same in 2 dets
         CALL TISET('GETHELEM2 ',ISUB)
         CALL SltCnd(nEl,nBasisMax,nBasis,NI,NJ,G1,nEl-iC,NMSH,FCK,TMat,NMAX,ALAT,UMat,Sum)
         GetHElement2=Sum
         IF(iC.EQ.0) GetHElement2%v=GetHElement2%v+ECore
         CALL TIHALT('GETHELEM2 ',ISUB)
         RETURN
      END

!  We always initialize from reals
       SUBROUTINE HDElemFromVal(h,v)
         TYPE(HDElement) h
         REAL*8  v
         h%v=v
         RETURN
      END
!.. Or HElements
       SUBROUTINE HDElemFromHElem(h,h2)
         TYPE(HDElement) h
         TYPE(HElement) h2
#ifdef __CMPLX
         IF(ABS(DIMAG(h2%v)).lt.1.D-10) THEN
            WRITE(6,*) "Conversion from complex to real:",h2%v
            STOP "Imaginary part non-zero"
         ENDIF
         h%v=DREAL(h2%v)
#else
         h%v=h2%v
#endif
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemAdd(h1,h2)
         TYPE(HDElement) h1,h2
         HDElemAdd%v=h1%v+h2%v
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemSub(h1,h2)
         TYPE(HDElement) h1,h2
         HDElemSub%v=h1%v-h2%v
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemNeg(h1)
         TYPE(HDElement) h1
         HDElemNeg%v=-h1%v
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemMul(h1,h2)
         TYPE(HDElement) h1,h2
         HDElemMul%v=h1%v*h2%v
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemDiv(h1,h2)
         TYPE(HDElement) h1,h2
         HDElemDiv%v=h1%v/h2%v
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemDivInt(h1,i2)
         TYPE(HDElement) h1
         INTEGER i2
         HDElemDivInt%v=h1%v/i2
         RETURN
      END

      TYPE(HDElement) FUNCTION HDElemExp(h1)
         TYPE(HDElement) h1
         HDElemExp%v=EXP(h1%v)
         RETURN
      END
      TYPE(HDElement) FUNCTION HDElemLog(h1)
         TYPE(HDElement) h1
         HDElemLog%v=LOG(h1%v)
         RETURN
      END
      REAL*8 FUNCTION HDElemSq(h1)
         TYPE(HDElement) h1
         HDElemSq=h1%v**2
         RETURN
      END

      LOGICAL FUNCTION HDElemAGT(h1,r2)
         TYPE(HDElement) h1
         REAL*8 r2
         HDElemAGT=ABS(h1%v).GT.r2
         RETURN
      END
      LOGICAL FUNCTION HDElemAGE(h1,r2)
         TYPE(HDElement) h1
         REAL*8 r2
         HDElemAGE=ABS(h1%v).GE.r2
         RETURN
      END
      REAL*8 FUNCTION HDElemDReal(h)
         TYPE(HDElement) h
         HDElemDReal=h%v
         RETURN
      END
END MODULE HElement

