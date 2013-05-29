! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      INTEGER FUNCTION ICMPDETS(D1,D2,NEL)
         IMPLICIT NONE
         INTEGER NEL,I
         INTEGER D1(1:NEL),D2(1:NEL),ICMPCSFS
         LOGICAL iscsf_old
         IF(iscsf_old(D1,NEL).OR.iscsf_old(D2,NEL)) THEN
            ICMPDETS=ICMPCSFS(D1,D2,NEL)
            RETURN
         ENDIF
         DO I=1,NEL
            IF(D1(I).LT.D2(I)) THEN
               ICMPDETS=-1
               RETURN
            ELSE IF(D1(I).GT.D2(I)) THEN
               ICMPDETS=1
               RETURN
            ENDIF
         ENDDO
         ICMPDETS=0
         RETURN
      END

      INTEGER FUNCTION IGETEXCITLEVEL(NI,NJ,NEL)
         IMPLICIT NONE
         INTEGER I,J,NEL,NI(NEL),NJ(NEL),IC
         LOGICAL iscsf_old
         IF(iscsf_old(NI,NEL).OR.iscsf_old(NJ,NEL)) THEN
            CALL GETCSFEXCITLEVEL(NI,NJ,NEL,IGETEXCITLEVEL)
C            IGETEXCITLEVEL=2
            RETURN
         ENDIF
         IC=0
         DO I=1,NEL
            DO J=1,NEL
               IF(NI(I).EQ.NJ(J)) THEN
                   IC=IC+1
                   EXIT
               ENDIF
            ENDDO
         ENDDO
         IGETEXCITLEVEL=NEL-IC
         RETURN
      END

!This routine is the same as IGetExcitLevel, but for two differences.
!First, this will only work for ordered lists of occupied orbitals
!Secondly, this will only find the number of excitation levels apart up
!to a maximum excitation level of MaxExcit. If the distance away is more than this
!it will simply return MaxExcit+1. i.e. if only want to know if it is a double
!excitation or less, then CALL IGetExcitLevel_2(NI,NJ,NEL,2). This will return
!0, 1, 2, or 3 if the distance is more than a double excitation.
      INTEGER FUNCTION IGetExcitLevel_2(NI,NJ,NEL,MaxExcit)
        IMPLICIT NONE
        INTEGER :: NEL,NI(NEl),NJ(NEl),MaxExcit,I,IC,j,MaxIC
        LOGICAL :: CommonOrb
        LOGICAL :: iscsf_old

        IF(iscsf_old(NI,NEL).OR.iscsf_old(NJ,NEL)) THEN
            CALL GETCSFEXCITLEVEL(NI,NJ,NEL,IGETEXCITLEVEL_2)
            RETURN
        ENDIF
        IF(MaxExcit.ge.NEl) THEN
            MaxIC=NEl-1
        ELSE
            MaxIC=MaxExcit
        ENDIF
        I=1   !Index for electron in NI
        IC=0  !Number of excitation levels between determinants
        DO WHILE((IC.le.MaxIC).and.(I.le.NEl))
            CommonOrb=.false.
            do j=1,NEL
                IF(NJ(j).eq.NI(I)) THEN
                    CommonOrb=.true.
                    EXIT
                ELSEIF(NJ(j).gt.NI(I)) THEN
                    EXIT
                ENDIF
            enddo
            IF(.not.CommonOrb) THEN
                IC=IC+1
            ENDIF
            I=I+1
        ENDDO
        IGetExcitLevel_2=IC
        RETURN
      END


C.. Very unnecessarily complex, but would be faster for large numbers of
C.. electrons.

      INTEGER FUNCTION IGETEXCITLEVEL_(NI,NJ,NEL)
         IMPLICIT NONE
         INTEGER I,J,NEL,NI(NEL),NJ(NEL),IC
         IC=0
         I=1
         J=1
C.. We only count differences from I to J
         DO WHILE(I.LE.NEL.AND.J.LE.NEL)
            DO WHILE(NI(I).LT.NJ(J).AND.I.LE.NEL)
               I=I+1
               IC=IC+1
            ENDDO
            DO WHILE(NI(I).GT.NJ(J).AND.I.LE.NEL.AND.J.LE.NEL)
               J=J+1
C               IC=IC+1
            ENDDO
            IF(NI(I).EQ.NJ(J)) THEN
               I=I+1
               J=J+1
            ENDIF            
         ENDDO
         IC=IC+(NEL+1-I)
         IGETEXCITLEVEL_=IC 
         RETURN
      END


! This function only finds the excitation level
! w.r.t. HF, by comparing the orbital energy to the chemical potential of
! the system. This can be made quicker by using an inverse BRR, which
! would only involve integer comparisons, but this is not exactly large
! speed increases...
      INTEGER FUNCTION iGetExcitLevelwrtHF(nI,NEl)
        USE IntegralsData , only : ChemPot
        USE SystemData , only : Arr
        IMPLICIT NONE
        INTEGER i,NEl,nI(NEl),IC
        
        IC=NEl
        do i=1,NEl
            IF(Arr(nI(i),2).lt.ChemPot) IC=IC-1
        enddo
        IGetExcitLevelwrtHF=IC

        RETURN

      END FUNCTION iGetExcitLevelwrtHF
      
      SUBROUTINE  CALCRHOIIMAX(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,ECORE)
         use constants, only: dp
         INTEGER IDET(NEL),I_P,NEL
         real(dp) :: Beta,RhoMax,ECore
         CALL CALCRHOIIMAXR(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,IDET,1,0,ECORE)

         RETURN
      END

      RECURSIVE SUBROUTINE CALCRHOIIMAXR(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,IDET,NTHISEL,ILASTBASIS,ECORE)
         use constants, only : dp
         INTEGER IDET(NEL),NTHISEL,I,ILASTBASIS
         real(dp) RH,RHOMAX,BETA,ECORE
         CHARACTER(*), PARAMETER :: t_r="CALCRHOIIMAXR"
         DO I=ILASTBASIS+1,NBASIS
            IDET(NTHISEL)=I
            call stop_all(t_r,"Inconsistent arguments")
            IF(NTHISEL.EQ.NEL) THEN
!               CALL CALCRHO2(IDET,IDET,BETA,I_P,NEL,NBASISMAX,
!     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,RH,NTAY,0,ECORE)
!               IF(RH.GT.RHOMAX) RHOMAX=RH
            ELSE
!               CALL CALCRHOIIMAXR(BETA,I_P,NEL,
!     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
!     &           RHOMAX,IDET,NTHISEL+1,I,ECORE)
            ENDIF
         ENDDO
         RETURN
      END
 
      SUBROUTINE  CALCRHOIIMAX2(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           ISTEPS,ECORE,RHOMAX)
         use constants, only : dp
         INTEGER NI(NEL)
         real(dp) RHOMAX
         INTEGER BRR(NBASIS)
         character(*), parameter :: t_r="CALCRHOIIMAX2"
         RHOMAX=0.0_dp
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         DO I=1,ISTEPS
             call stop_all(t_r,"Call to GenNewDet2 removed due to       &
     &          incorrect number of arguments")
!            CALL GENNEWDET2(NI,NI,NEL,NBASIS,NBASISMAX,G1,ISEED)
!            CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,
!     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,RH,NTAY,0,ECORE)
!            IF(RH.GT.RHOMAX) RHOMAX=RH
         ENDDO
         RETURN
      END


C.. GETHELEMENT
C.. Get matrix element of the hamiltonian
      FUNCTION GETHELEMENT(II,IJ,HAMIL,LAB,NROW,NDET)
      use constants, only: dp
      IMPLICIT NONE
      real(dp) HAMIL(*),GETHELEMENT
      INTEGER LAB(*)
      INTEGER NDET,NROW(NDET),IJ,II,I,J,INDXROW,IMAX,K
C.. We only have half of H, so if J<I, return the symmetrical (J,I) element
C.. Or if we have the whole H, it's quicker to look closer to its beginning
      IF(IJ.LT.II) THEN
         I=IJ
         J=II
      ELSE
         I=II
         J=IJ
      ENDIF
      GETHELEMENT=0.0_dp
      INDXROW=1
C.. Find the Ith row
      DO K=1,I-1
         INDXROW=INDXROW+NROW(K)
      ENDDO
      IMAX=INDXROW+NROW(I)-1
      DO K=INDXROW,IMAX
         IF(LAB(K).GT.J) RETURN
         IF(LAB(K).EQ.J) THEN
            GETHELEMENT=HAMIL(K)
            RETURN
         ENDIF
      ENDDO
      RETURN
      END

C.. GETHROW
C.. Get the index within HAMIL of a row II of the H matrix
      INTEGER FUNCTION IGETHROW(II,NROW,NDET)
      IMPLICIT NONE
      INTEGER NDET,NROW(NDET),K,II
      IGETHROW=1
C.. Find the Ith row
      DO K=1,II-1
         IGETHROW=IGETHROW+NROW(K)
      ENDDO
      RETURN
      END

