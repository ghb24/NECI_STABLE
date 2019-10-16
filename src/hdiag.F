

!Diagonalize a compressed hamiltonian (in HAMIL), returning eigenvalues in W and eigenvectors in CK.
       SUBROUTINE HDIAG_neci(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,
     &                     NBLOCKSTARTS,NBLOCKS)
         use constants, only: dp,int32!,n_int
         use SystemData, only: BasisFN
         use global_utilities
         use HElem, only: HElement_t_size
!         use DetCalcData, only: NMRKS
!         use DetBitOps, only: count_open_orbs,EncodeBitDet,DetBitEQ
!         use DetBitOps, only: FindBitExcitLevel
!         use Bit_Reps
!         use HPHFRandExcitMod, only: FindExcitBitDetSym
         IMPLICIT NONE
!         INTEGER(n_int) :: BitRep(0:NIfTot),BitRepSym(0:NIfTot)
!         INTEGER(n_int) :: TempBitRep(0:NIfTot),BitRepHF(0:NIfTot)
         INTEGER NDET
         INTEGER NROW(NDET),LAB(*)
         HElement_t(dp) HAMIL(*),CK(NDET,NDET),WORK2(3*NDET),WORK
!         HElement_t(dp) CopyCK(NDET,NDET)
         real(dp) W(NDET)
         INTEGER I,J,IND,INDZ,NBS!,nopen,k,Excit
         INTEGER(int32) INFO
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS+1)
         type(timer), save :: proc_timer
         real(dp) GSEN
         character(len=*), parameter :: t_r="HDIAG_neci"
         proc_timer%timer_name='HDIAG     '
         call set_timer(proc_timer)
         GSEN=1.D100
         DO I=1,NDET
            DO J=1,NDET
               CK(I,J)=0.0_dp
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
         IF(HElement_t_size.ne.1) THEN
             !ensure hermitian
             do i=1,ndet
                 do j=1,i-1
#ifdef __CMPLX
                     !To avoid compiler warnings
                     CK(i,j)=CONJG(CMPLX(CK(j,i),kind=kind(0.0_dp)))
#endif
                 enddo
             enddo
         ENDIF

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
            IF(HElement_t_size.EQ.1) THEN
#ifndef __CMPLX
               CALL DSYEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK2,3*NDET,INFO)
#else
               call stop_all(t_r,"This block for real calcs only")
#endif
            ELSE
#ifdef __CMPLX
                 CALL ZHEEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK,4*NDET,WORK2,INFO)
#else
                 call stop_all(t_r,"This block for complex calcs only")
#endif
            ENDIF
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               call stop_all(t_r, "DSYEV error")
            ENDIF
            IF(W(NBS).LT.GSEN) GSEN=W(NBS)
! This is commented out develomental code to label the determinants
! by their HPHF partner and number of open shells.
!
!            do j=1,ndet
!                call EncodeBitDet(NMRKS(:,j),BitRep)
!                if(j.eq.1) then
!                    BitRepHF=BitRep
!                endif
!                nopen = count_open_orbs(BitRep)
!                Excit = FindBitExcitLevel(BitRepHF,BitRep,NEl)
!                call FindExcitBitDetSym(BitRep,BitRepSym)
!                do k=1,ndet
!                    call EncodeBitDet(NMRKS(:,k),TempBitRep)
!                    if(DetBitEQ(TempBitRep,BitRepSym,NIfDBO)) then
!!                        write(6,*) "Found same: ",j,k,BitRep(0:NIfDBO), &
!!     &                      TempBitRep(0:NIfDBO)
!                        exit
!                    endif
!                enddo
!                if(k.eq.ndet+1) call stop_all(this_routine, 'Error in finding spin-coup dets')
!                write(69,"(I4,4G20.10,3I5)") j,CK(j,1:4),Excit,nopen,k
!            enddo
!
         ENDDO

!Developmental and research routine
!         CALL ExcitChangeTest(CK,W,CopyCK,NDET)

C.. CK now contains the eigenvectors, and W the eigenvalues
         WRITE(6,"(A,F19.11,I4)") "GROUND E=",GSEN
         call halt_timer(proc_timer)
         RETURN
      END

!This routine for each determinant, calculates each excitation, and then its contribution to the change in
!the coefficient of the desired eigenvector. The average for this is obviously zero.
!This is called from HDIAG, but will generally be commented out, since it is just for interest.
      SUBROUTINE ExcitChangeTest(Vecs,Vals,HMAT,NDET)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NDET,i,j,NExcit,bin,nbins
         HElement_t(dp) Vecs(NDET,NDET),HMAT(NDET,NDET),KMAT(NDET,NDET)
         real(dp) Vals(NDET),S(NDET,NDET+1),SSUM
         real(dp) MinS,MaxS,RangeBins
         INTEGER , ALLOCATABLE :: SHist(:,:)
         real(dp) , ALLOCATABLE :: RangeHist(:)
         character(*), parameter :: t_r = 'ExcitChangeTest'

!Vals(1) contains ground state
!Check that eigenvectors are correct by multiplying matrix by first eigenvector
!         TrialEigenvec=0.0_dp
!         Norm=(0.0_dp)
!         do i=1,NDET
!             do j=1,NDET
!                 TrialEigenvec(i)=TrialEigenvec(i)+HMAT(i,j)*Vecs(j,1)
!             enddo
!             Norm=Norm+(Vecs(i,1)*Vecs(i,1))
!         enddo
!         WRITE(6,*) "Norm for first eigenvector is: ", Norm
!         do i=1,NDET
!             TrialEigenvec(i)=TrialEigenvec(i)/(Vals(1))
!             WRITE(6,*) Vals(i)
!             WRITE(6,*) TrialEigenvec(i),Vecs(i,1)
!         enddo

!Create the k matrix...
         KMAT=(0.0_dp)
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
         S=0.0_dp
         MinS=100.0
         MaxS=-100.0
         WRITE(11,*) "# Det   SijSum     MinS    MaxS"
         do i=1,NDET

!Count excitations attached to it...
             NExcit=0
             do j=1,NDET
!We only want excitations, so remove the diagonal contribution
                 IF(i.eq.j) CYCLE
                 IF(abs(KMAT(i,j)).gt.0.0_dp) NExcit=NExcit+1
             enddo

!             IF(NExcit.lt.1) call stop_all(this_routine, 'ERROR IN COUNTING EXCITS')

             do j=1,NDET
                 IF(i.eq.j) CYCLE
                 IF(abs(KMAT(i,j)).gt.0.0_dp) THEN
!For each excitation of each determinant i, find its contribution to the change of expansion coefficient
                     S(i,j)=((KMAT(i,i)/(NExcit+0.0_dp))*Vecs(i,1))+
     &                        ((KMAT(i,j))*Vecs(j,1))
                     IF(S(i,j).lt.MinS) MinS=S(i,j)
                     IF(S(i,j).gt.MaxS) MaxS=S(i,j)
                 ENDIF
             enddo

!Sum all the changes of the expansion coefficient for a given i due to the excitations j
             SSUM=0.0_dp
             do j=1,NDET
                 SSUM=SSUM+S(i,j)
             enddo
!Store this value in the final column of the S array
             WRITE(6,"(I4,3G25.16)") i,SSUM,MinS,MaxS
             WRITE(11,"(I4,3G25.16)") i,SSUM,MinS,MaxS
             S(i,NDET+1)=SSUM

         enddo

         RangeBins=(MaxS-MinS)/(NBins+0.0_dp)

         ALLOCATE(SHist(NDET,NBins))
         SHist(1:nDet,1:NBins)=0
         ALLOCATE(RangeHist(NBins+1))
         RangeHist=0.0_dp

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
                 IF(abs(KMAT(i,j)).gt.0.0_dp) THEN
!j is excitation of i - find which bin the overlap goes into
                     Bin=1
                     do while(S(i,j).ge.RangeHist(Bin+1))
                         Bin=Bin+1
                     enddo
                     IF(Bin.gt.NBins)
     &                  call stop_all(t_r,'Error in creating histogram')
                     SHist(i,Bin)=SHist(i,Bin)+1

                 ENDIF
             enddo

         enddo

         do i=1,NBins
            WRITE(9,"(G20.12)",advance='no')
     &                     (MinS+((i-1.0_dp)*RangeBins))+(RangeBins/2.0)
            do j=1,NDET
                WRITE(9,"(I8)",advance='no') SHist(j,i)
            enddo
            WRITE(9,*) ""
        enddo

        DEALLOCATE(SHist)
        DEALLOCATE(RangeHist)

      END SUBROUTINE ExcitChangeTest
