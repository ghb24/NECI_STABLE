!Diagonalize a compressed hamiltonian (in HAMIL), returning eigenvalues in W and eigenvectors in CK.
       SUBROUTINE HDIAG(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,LENHAMIL,
     &                     NBLOCKSTARTS,NBLOCKS,BLOCKSYM)
         USE HElem
         use System, only: BasisFN
         IMPLICIT NONE
         INTEGER NDET,LENHAMIL
         INTEGER NROW(NDET),LAB(*)
         TYPE(HElement)  HAMIL(*),CK(NDET,NDET),WORK2,WORK
!         TYPE(HElement) CopyCK(NDET,NDET)
         REAL*8 W(NDET)
         INTEGER I,J,IND,INDZ,INFO,LDA,NBS
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS),ISUB
         TYPE(BasisFN) BLOCKSYM(NBLOCKS)
         REAL*8 GSEN
         CALL TISET('HDIAG     ',ISUB)
         GSEN=1.D100
         DO I=1,NDET
            DO J=1,NDET
               CK(I,J)=0.D0
            ENDDO
         ENDDO
C.. Now we fill the RIJ array
         IND=1
         INDZ=1
         DO I=1,NDET
            INDZ=INDZ+NROW(I)
            DO WHILE (IND.LT.INDZ)
               CK(I,LAB(IND))=HAMIL(IND)
               IND=IND+1
            ENDDO
         ENDDO
!****************************
!Needs to be uncommented for ExcitChangeTest to work
!         DO I=1,NDET
!             DO J=I,NDET
!                 CopyCK(I,J)=CK(I,J)
!                 CopyCK(J,I)=CK(I,J)
!             ENDDO
!         ENDDO
!         DO I=1,NDET
!             DO J=1,NDET
!                 WRITE(6,"(F20.14)",advance='no') CopyCK(I,J) 
!             ENDDO
!             WRITE(6,*) ""
!         ENDDO
!***************************

C.. Diagonalize
         DO I=1,NBLOCKS
            NBS=NBLOCKSTARTS(I)
            IF(HElementSize.EQ.1) THEN
               CALL DSYEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK2,3*NDET,INFO)
            ELSE
                 CALL ZHEEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK,4*NDET,WORK2,INFO)
            ENDIF
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
            IF(W(NBS).LT.GSEN) GSEN=W(NBS)
         ENDDO

!Developmental and research routine
!         CALL ExcitChangeTest(CK,W,CopyCK,NDET)

C.. CK now contains the eigenvectors, and W the eigenvalues
         WRITE(6,"(A,F19.11,I4)") "GROUND E=",GSEN
         CALL TIHALT('HDIAG     ',ISUB)
         RETURN
      END

!This routine for each determinant, calculates each excitation, and then its contribution to the change in
!the coefficient of the desired eigenvector. The average for this is obviously zero.
!This is called from HDIAG, but will generally be commented out, since it is just for interest.
      SUBROUTINE ExcitChangeTest(Vecs,Vals,HMAT,NDET)
         USE HElem
         IMPLICIT NONE
         INTEGER NDET,i,j,NExcit,INFO,bin,nbins
         TYPE(HElement) Vecs(NDET,NDET),HMAT(NDET,NDET),KMAT(NDET,NDET)
         TYPE(HElement) TrialEigenvec(NDET),Norm
         REAL*8 Vals(NDET),KVals(NDET),Work(3*NDET),S(NDET,NDET+1),SSUM
         REAL*8 MinS,MaxS,RangeBins
         INTEGER , ALLOCATABLE :: SHist(:,:)
         REAL*8 , ALLOCATABLE :: RangeHist(:)

!Vals(1) contains ground state
!Check that eigenvectors are correct by multiplying matrix by first eigenvector
!         CALL AZZERO(TrialEigenvec,NDET*HElementSize)
!         Norm=HElement(0.D0)
!         do i=1,NDET
!             do j=1,NDET
!                 TrialEigenvec(i)=TrialEigenvec(i)+HMAT(i,j)*Vecs(j,1)
!             enddo
!             Norm=Norm+(Vecs(i,1)*Vecs(i,1))
!         enddo
!         WRITE(6,*) "Norm for first eigenvector is: ", Norm%v
!         do i=1,NDET
!             TrialEigenvec(i)=TrialEigenvec(i)/HElement(Vals(1))
!             WRITE(6,*) Vals(i)
!             WRITE(6,*) TrialEigenvec(i),Vecs(i,1)
!         enddo

!Create the k matrix...
         CALL AZZERO(KMAT,NDET*HElementSize)
         DO i=1,NDET
             KMAT(i,i)=Vals(1)
         ENDDO
         DO i=1,NDET
             DO j=1,NDET
                 KMAT(i,j)=KMAT(i,j)-HMAT(i,j)
             ENDDO
         ENDDO

!         DO I=1,NDET
!             DO J=1,NDET
!                 WRITE(6,"(F20.14)",advance='no') KMAT(I,J) 
!             ENDDO
!             WRITE(6,*) ""
!         ENDDO

!Diagonalise to test eigenvectors are the same
!         CALL DSYEV('V','U',NDET,KMAT,
!     &         NDET,KVals,Work,3*NDET,INFO)

!         DO i=1,NDET
!             WRITE(6,*) Vals(i),KVals(i),Vals(1)-Vals(i)
!         ENDDO

!We now want eigenvector corresponding to largest eigenvalue of KMAT or smallest of HMAT
!         DO i=1,NDET
!             WRITE(6,*) Vecs(i,1), KMAT(i,NDET)
!         ENDDO
         

!Run through determinants in graph
         NBins=100
         CALL AZZERO(S,NDET*(NDET+1))
         MinS=100.0
         MaxS=-100.0
         WRITE(11,*) "# Det   SijSum     MinS    MaxS"
         do i=1,NDET

!Count excitations attached to it...
             NExcit=0
             do j=1,NDET
!We only want excitations, so remove the diagonal contribution
                 IF(i.eq.j) CYCLE
                 IF(KMAT(i,j).agt.0.D0) NExcit=NExcit+1
             enddo

!             IF(NExcit.lt.1) STOP 'ERROR IN COUNTING EXCITS'

             do j=1,NDET
                 IF(i.eq.j) CYCLE
                 IF(KMAT(i,j).agt.0.D0) THEN
!For each excitation of each determinant i, find its contribution to the change of expansion coefficient
                     S(i,j)=((KMAT(i,i)%v/(NExcit+0.D0))*Vecs(i,1)%v)+
     &                        ((KMAT(i,j)%v)*Vecs(j,1)%v)
                     IF(S(i,j).lt.MinS) MinS=S(i,j)
                     IF(S(i,j).gt.MaxS) MaxS=S(i,j)
                 ENDIF
             enddo

!Sum all the changes of the expansion coefficient for a given i due to the excitations j
             SSUM=0.D0
             do j=1,NDET
                 SSUM=SSUM+S(i,j)
             enddo
!Store this value in the final column of the S array
             WRITE(6,"(I4,3G25.16)") i,SSUM,MinS,MaxS
             WRITE(11,"(I4,3G25.16)") i,SSUM,MinS,MaxS
             S(i,NDET+1)=SSUM

         enddo

         RangeBins=(MaxS-MinS)/(NBins+0.D0)

         ALLOCATE(SHist(NDET,NBins))
         CALL IAZZERO(SHist,NBins*NDET)
         ALLOCATE(RangeHist(NBins+1))
         CALL AZZERO(RangeHist,NBins+1)

         RangeHist(1)=MinS
         do i=2,NBins
             RangeHist(i)=RangeHist(i-1)+RangeBins
         enddo
!Make sure final bin catches them all
         RangeHist(NBins+1)=MaxS+10

         do i=1,NDET
!Run through determinants again, and determinant non-zero matrix elements

             do j=1,NDET
                 IF(i.eq.j) CYCLE
                 IF(KMAT(i,j).agt.0.D0) THEN
!j is excitation of i - find which bin the overlap goes into
                     Bin=1
                     do while(S(i,j).ge.RangeHist(Bin+1))
                         Bin=Bin+1
                     enddo
                     IF(Bin.gt.NBins) STOP 'Error in creating histogram'
                     SHist(i,Bin)=SHist(i,Bin)+1

                 ENDIF
             enddo

         enddo

         do i=1,NBins
            WRITE(9,"(G20.12)",advance='no') (MinS+((i-1.D0)*RangeBins))
     &                                    +(RangeBins/2.0)
            do j=1,NDET
                WRITE(9,"(I8)",advance='no') SHist(j,i)
            enddo
            WRITE(9,*) ""
        enddo

        DEALLOCATE(SHist)
        DEALLOCATE(RangeHist)

      END SUBROUTINE ExcitChangeTest

 
      SUBROUTINE HDIAG_NH(NDET,NBLOCKSTARTS,NBLOCKS,
     &            NEL,NMRKS,NBASISMAX,NBASIS,G1,NMSH,BRR,
     &            FCK,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,
     &            BETA,I_P,ILOGGING,IFDET,ARR,BLOCKSYM)
         use System, only: BasisFN
         IMPLICIT NONE
         INTEGER NDET,NEL,NBLOCKS,NBLOCKSTARTS(NBLOCKS),NMRKS(NEL,NDET)
         INTEGER NBASISMAX(5,2),NBASIS
         TYPE(BasisFN) G1(*)
         INTEGER NMSH,NMAX,I_P,ILOGGING,BRR(*),INFO
         REAL*8  ALAT(3),UMAT(*),ICMAX(*),GC(*),BETA
         LOGICAL TMC,TWARN
         REAL*8  FCK(*),ECORE,GETHELEMENT2,DNORM
         REAL*8  CALCDLWDB,CALCMCEN,ARR(NBASIS)
         
         INTEGER I,NBLKMAX,NBLK,J,K,NBS,III,IFDET,L
         REAL*8 CK(*),W(*),WORK(*),EXEN,GSEN,FLRI,FLSI,GSEN2
         POINTER (IP_CK,CK),(IP_W,W),(IP_WORK,WORK)
         REAL*8 ELOWEST,GOFE(NBASIS),EBIAS
         INTEGER GELOWEST,ISUB,IDEG
         TYPE(BasisFN) BLOCKSYM(NBLOCKS)
         CALL TISET('HDIAG_NH  ',ISUB)
C.. Find the largest block
         NBLKMAX=0
         DO I=1,NBLOCKS
            IF(NBLOCKSTARTS(I+1)-NBLOCKSTARTS(I).GT.NBLKMAX)
     &         NBLKMAX=NBLOCKSTARTS(I+1)-NBLOCKSTARTS(I)
         ENDDO
         WRITE(6,*) "Maximum Block Size:",NBLKMAX
C.. Now allocate memory
         CALL MEMORY(IP_CK,NBLKMAX*NBLKMAX,"CK")
         CALL MEMORY(IP_W,NBLKMAX,"W")
         CALL MEMORY(IP_WORK,3*NBLKMAX,"WORK")
C.. Go through the blocks forming the hamil matrix
         EXEN=0.D0
         DNORM=0.D0
         EBIAS=0.D0
         DO J=1,NBASIS
            GOFE(J)=0.D0
         ENDDO
         OPEN(14,FILE='RHOPIIex',STATUS='UNKNOWN')
         OPEN(15,FILE='ENERGIES',STATUS='UNKNOWN')
         ELOWEST=1.D200
         GELOWEST=0
         DO I=0,NBLOCKS-1
            NBS=NBLOCKSTARTS(I+1)
            NBLK=NBLOCKSTARTS(I+2)-NBS
            DO J=0,NBLK-1
               DO K=J,NBLK-1
                  CK(J*NBLK+K+1)=GETHELEMENT2(NMRKS(1,NBS+J),
     &                      NMRKS(1,NBS+K),NEL,NBASISMAX,
     &   G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,-1,ECORE)
               ENDDO
            ENDDO
C.. Diagonalize
            CALL DSYEV('V','L',NBLK,CK,NBLK,W,WORK,3*NBLKMAX,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
C.. Now work out exact wi etc.
            CALL GETSYMDEGEN(BLOCKSYM(I+1),NBASISMAX,IDEG)
            DO III=1,NBLK
               WRITE(15,"(2I4,I7,$)") I+1,III,NBS+III-1
!               CALL WRITEDET(15,NMRKS(1,NBS+III-1),NEL,.FALSE.)
               WRITE(15,"(F19.11,I5)") W(III),IDEG

               IF(ABS(W(III)/ELOWEST-1.D0).LT.1.D-6) THEN
                  GELOWEST=GELOWEST+IDEG
               ELSEIF(W(III).LT.ELOWEST) THEN
C.. if we're going to change ELOWEST, and there was a previous ELOWEST,
C.. we need to multiply DNORM and EXEN by exp(-BETA(ENEW-EOLD))
                     EXEN=EXEN*EXP(-BETA*(EBIAS-W(III)))
                     DNORM=DNORM*EXP(-BETA*(EBIAS-W(III)))
                  ELOWEST=W(III)
                  GELOWEST=IDEG
                  EBIAS=W(III)
               ENDIF
C
C               DO J=1,NBLK
C                  WRITE(15,*) CK((III-1)*NBLK+J)
C               ENDDO
               CALL CALCRHOPII(III,NBLK,NBLK,CK,W,BETA,I_P,ILOGGING,
     &               ECORE,FLRI,FLSI,TWARN)
               WRITE(14,"(I12,A,$)") NBS+III-1,"("
               DO L=1,NEL
                  WRITE(14,"(I3,A,$)") NMRKS(L,NBS+III-1),","
               ENDDO
               GSEN=CALCDLWDB(III,NBLK,NBLK,CK,W,BETA,ECORE)
               WRITE(14,"(A,4G25.16,I5)") ") ",
     &            EXP(FLSI+I_P*FLRI),FLRI*I_P,FLSI,GSEN,IDEG
               EXEN=EXEN+(W(III)*IDEG)*EXP(-(W(III)-EBIAS)*BETA)
               DNORM=DNORM+IDEG*EXP(-(W(III)-EBIAS)*BETA)
C.. Find the lowest energy, and write out g(E)
C               IF(W(III).LT.ELOWEST) THEN
C                  ELOWEST=W(III)
                  DO J=1,NBLK
                     DO K=1,NEL
                        L=NMRKS(K,NBS+J-1)
                       GOFE(L)=GOFE(L)
     &                   +IDEG*(CK((III-1)*NBLK+J)**2)*EXP(-W(III)*BETA)
C                        IF(L.EQ.9.AND.
C     &                  ABS(CK((III-1)*NBLK+J)).GT.1.D-5) THEN
C         WRITE(6,*) NBS+III-1,J,CK((III-1)*NBLK+J),W(III),
C     &         (CK((III-1)*NBLK+J)**2)*EXP(-W(III)*BETA)
C                        ENDIF
                     ENDDO
                  ENDDO
C               ENDIF
            ENDDO
            IF(IFDET.GE.NBS.AND.IFDET.LT.NBLOCKSTARTS(I+2)) THEN
               GSEN2=CALCDLWDB(IFDET-NBS+1,NBLK,NBLK,CK,W,BETA,ECORE)
            ENDIF
         ENDDO
         CLOSE(14)
         CLOSE(15)
C        WRITE(6,*) "NORM: ",DNORM
         WRITE(6,"(A,F19.11)") "EXACT E(BETA)=",EXEN/DNORM
         WRITE(6,"(A,F19.11)") "EXACT DLWDB(D0)=",GSEN2
         WRITE(6,"(A,F19.11,I4)") "GROUND E, g(E)=",ELOWEST,GELOWEST
         OPEN(16,FILE='GOFE',STATUS='UNKNOWN')
         DO J=1,NBASIS
            WRITE(16,"(6I4,2G25.16)") BRR(J),(G1(BRR(J))%k(K),K=1,5),
     &         ARR(J),GOFE(BRR(J))/DNORM
         ENDDO
         CLOSE(16)
         CALL FREEM(IP_WORK)
         CALL FREEM(IP_CK)
         CALL FREEM(IP_W)
         CALL TIHALT('HDIAG_NH  ',ISUB)
         RETURN
      END 
