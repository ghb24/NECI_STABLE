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
#ifdef __CMPLX
         HElemDConjg%v=DCONJG(h%v)
#else
         HElemDConjg%v=h%v
#endif
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
      TYPE(HElement) FUNCTION HElemPow(h,r)
         TYPE(HElement) h
         REAL*8 r
         HElemPow%v=h%v**r
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
      TYPE(HDElement) FUNCTION HElemABS(h1)
         TYPE(HElement) h1
         HElemABS%v=ABS(h1%v)
         RETURN
      END
      REAL*8 FUNCTION HElemSq(h1)
         TYPE(HElement) h1
#ifdef __CMPLX
         HElemSq=h1%v*dconjg(h1%v)
#else
         HElemSq=h1%v**2
#endif
         RETURN
      END
      REAL*8 FUNCTION HElemDReal(h1)
         TYPE(HElement) h1
#ifdef __CMPLX
         HElemDReal=DREAL(h1%v)
#else
         HElemDReal=h1%v
#endif
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

      TYPE(HElement) FUNCTION GetHElement2(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,iC2,ECore)
         INTEGER NMSH,NMAX
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*)
         INCLUDE 'basis.inc'
         INCLUDE 'uhfdet.inc'
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(*)
         TYPE(HElement) UMat(*)
         INTEGER I,nEl,NI(nEl),NJ(nEl),iC,nBasisMax(5,2),iC2
         REAL*8 ECore
         TYPE(HElement) Sum,Sum2
         INTEGER IGETEXCITLEVEL
         LOGICAL ISCSF
         INTEGER ISUB
         IF(ISCSF(NI,NEL).OR.ISCSF(NJ,NEL)) THEN
            CALL CSFGETHELEMENT(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,ECore,Sum2)
            GETHELEMENT2=SUM2
            RETURN
         ENDIF
         IF(tStoreAsExcitations.AND.nI(1).eq.-1.and.nJ(1).eq.-1) then
            if(ic2.ne.2) stop 'tStoreAsExcitations in GetHElement2 requires ic=2 (doubles).'
            Call SCR2Excit(nBasisMax,nJ,G1,nBasis,UMat,Alat,nBasisMax(2,3),Sum)
            GetHElement2=Sum
            RETURN
         endif
         IC=IC2
         GetHElement2%v=0.D0
         IF(IC.LT.0) IC=IGETEXCITLEVEL(NI,NJ,NEL)
!.. if we differ by more than 2 spin orbital, then the hamiltonian element is 0         
         IF(IC.GT.2) RETURN
!.. SLTCND has IC is # electrons the same in 2 dets
         CALL TISETL('GETHELEM2 ',ISUB,60)
         CALL SltCnd(nEl,nBasisMax,nBasis,NI,NJ,G1,nEl-iC,NMSH,FCK,NMAX,ALAT,UMat,Sum)
         GetHElement2=Sum
         IF(iC.EQ.0) GetHElement2%v=GetHElement2%v+ECore
!         CALL WRITEDET(6,NI,NEL,.FALSE.)
!         CALL WRITEDET(6,NJ,NEL,.FALSE.)
!         WRITE(6,*) GetHElement2
         CALL TIHALTL('GETHELEM2 ',ISUB,60)
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

!  Get a matrix element of the unperturbed Hamiltonian.  This is just the sum of the Hartree-Fock eigenvalues
      subroutine GetH0Element(nI,nEl,Arr,nBasis,ECore,hEl)
         use HElement
         implicit none
         integer nI(nEl),nEl,nBasis
         type(HElement) hEl
         real*8 Arr(nBasis,2),ECore
         integer i
         INCLUDE 'uhfdet.inc'
         if(tStoreAsExcitations.and.nI(1).eq.-1) then
!The excitation storage starts with -1.  The next number is the excitation level,L .  
!Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            hEl=0.d0
            do i=4,nI(2)+4-1
               hEl=hEl-HElement(Arr(nI(i),2))
            enddo
            do i=i,i+nI(2)-1
               hEl=hEl+HElement(Arr(nI(i),2))
            enddo
         else
            hEl=ECore
            do i=1,nEl
               hEl=hEl+HElement(Arr(nI(i),2))
            enddo
         endif
!         call writedet(77,nI,nel,.false.)
!         write(77,*) "H0",hEl
!         call flush(77)
      end

