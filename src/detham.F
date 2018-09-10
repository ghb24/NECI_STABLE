      SUBROUTINE DETHAM(NDET,NEL,NMRKS, HAMIL,LAB,NROW,TCOUNT, 
     &  ICMAX,GC,TMC)
      use global_utilities
      use constants, only: dp, n_int, eps
      Use Determinants, only: get_helement
      use SystemData , only : BasisFN, t_lattice_model, tGUGA
      use CalcData , only : TStar
      use bit_rep_data, only: NIfTot
      use SystemData, only: tHPHF
      use hphf_integrals, only: hphf_diag_helement
      use hphf_integrals, only: hphf_off_diag_helement
      use DetBitOps, only: EncodeBitDet

#ifndef __CMPLX
      use guga_matrixElements, only: calcDiagMatEleGuga_nI
      use guga_excitations, only: calc_guga_matrix_element, 
     &                            init_csf_information
      use guga_data, only: excitationInformation
#endif
      use lattice_mod, only: get_helement_lattice
      IMPLICIT NONE
      HElement_t(dp) HAMIL(*)
      INTEGER NDET,NEL
      INTEGER LAB(*)
      INTEGER NMRKS(NEL,*)
      INTEGER NROW(NDET),GC
      INTEGER ICMAX,KI,IBEG,IBEGJ,KJ,IMAX,IDAMAX
      LOGICAL TCOUNT,TMC
      INTEGER IGETEXCITLEVEL_2,STEP,IMAXJ
      integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
      HElement_t(dp) totSUM
      type(timer), save :: proc_timer
#ifndef __CMPLX 
      type(excitationInformation) :: excitInfo
#endif
      integer :: ic=0
!      LOGICAL TSTAR
C ==-------------------------------------------------------------------==
      proc_timer%timer_name='    DETHAM'
      call set_timer(proc_timer)
C ==-------------------------------------------------------------------==
C..Global counter
      GC=0
      NROW(1:NDET)=0
      IBEG=0
C..   Now we need to match up any two determinants
      DO KI=1,NDET
         call EncodeBitDet(NMRKS(:,KI), ilutI)
         IF(mod(KI,1000).eq.0) WRITE(6,*) KI
         IF(KI.EQ.1) THEN
           IBEG=0
         ELSE
           IBEG=IBEG+NROW(KI-1)
         ENDIF
         IF(TMC) THEN 
           IBEGJ=1
           STEP=1
           IMAXJ=NDET
         ELSE
           IF(TSTAR) THEN
               IF(KI.EQ.1) THEN
                  IBEGJ=1
                  IMAXJ=NDET
                  STEP=1
               ELSE
                  IBEGJ=KI
                  IMAXJ=KI
                  STEP=1
               ENDIF
           ELSE
               IBEGJ=KI
               IMAXJ=NDET
               STEP=1
           ENDIF
         ENDIF
         IF(STEP.EQ.0) STEP=1
         ! try to use the new guga matrix element calc here.. for this 
         ! misuse the current_stepvector quantities to init this only
         ! once for nI 
#ifndef __CMPLX 
         if (tGUGA) call init_csf_information(ilutI)
#endif
         DO KJ=IBEGJ,IMAXJ,STEP
            call EncodeBitDet(NMRKS(:,KJ), ilutJ)
            if (tHPHF) then
              if (KI == KJ) then
                totsum = hphf_diag_helement(NMRKS(:,KI), ilutI)
              else
                totsum = hphf_off_diag_helement(NMRKS(:,KI), 
     &                            NMRKS(:,KJ), ilutI, ilutJ)
              end if

#ifndef __CMPLX
            else if (tGUGA) then
              if (KI == KJ) then
                  totsum = calcDiagMatEleGuga_nI(NMRKS(:,KI))
              else
!                   totsum = calc_off_diag_guga_gen(ilutI, ilutJ)
                  call calc_guga_matrix_element(ilutI,ilutJ,excitInfo,
     &              totsum, .true., 1)
              end if
#endif

            else if (t_lattice_model) then 
              if (KI == KJ) then
                totsum=get_helement_lattice(NMRKS(:,KI),NMRKS(:,KJ),ic)
              else
                totsum=get_helement_lattice(NMRKS(:,KI),NMRKS(:,KJ))
              end if
            else 
              if (KI == KJ) then
                totsum = get_helement(NMRKS(:,KI), NMRKS(:,KJ), 0)
              else
                totsum = get_helement(NMRKS(:,KI), NMRKS(:,KJ),
     &                                 ilutI, ilutJ)
              end if
            end if
C            CALL SLTCND(NEL,NBASISMAX,NHG,
C     &     FCK,NMAX,ALAT,UMAT,TOTSUM)
            IF(ABS(TOTSUM).LT.1.0e-10_dp) TOTSUM=0.0_dp
            IF(abs(TOTSUM)>EPS .OR. (KI.EQ.KJ)) THEN
               GC=GC+1
C..   Stores the number of non-zero elements in each row
               NROW(KI)=NROW(KI)+1
               IF(.NOT.TCOUNT) THEN
                 LAB(IBEG+NROW(KI))=KJ
                 HAMIL(IBEG+NROW(KI))=TOTSUM
               ENDIF
            ENDIF
 500        CONTINUE
         ENDDO
      ENDDO
C..No. of columns
      IF(TCOUNT) THEN
!        IMAX=MAXLOC(NROW)
        IMAX=IDAMAX(NDET,real(NROW,dp),1)
        ICMAX=NROW(IMAX)
        WRITE(6,*) ' MAXIMUM WIDTH OF HAMIL : ' , ICMAX
        WRITE(6,*) ' TOTAL NUMBER OF NON-ZERO ELEMENTS : ' , GC
      ENDIF
C      IF(.NOT.TCOUNT) THEN
C        WRITE(6,*)
C        WRITE(6,*) ' SECTION OF MAP OF HAMILTONIAN : '
C        WRITE(6,*)
C        WRITE(6,*) NDET,ICMAX
C        DO I=1,NDET
C          WRITE(6,*) (LAB(J,I),J=1,ICMAX)
C          WRITE(6,*) I,(HAMIL(J,I),J=1,ICMAX)
C        ENDDO
C      ENDIF
C ==-------------------------------------------------------------------==
      call halt_timer(proc_timer)
C ==-------------------------------------------------------------------==
      RETURN
      END
C ==---------------------------------------------------------------==
      SUBROUTINE CHK_MS(ND1,ND2,IFLG,G1)
      use SystemData, only: BasisFN, nel
      use constants, only: dp
      IMPLICIT NONE
C..
      TYPE(BASISFN) G1(*)
      integer, intent(in) :: iFlg
      integer, intent(in) :: ND1(NEL), ND2(NEL)
      integer :: i, MS1, MS2
      character(*), parameter :: t_r = 'CHK_MS'
C..
      MS1=0
      MS2=0
      DO I=1,NEL
        MS1=MS1+G1(ND1(I))%MS
        MS2=MS2+G1(ND2(I))%MS
      ENDDO
      MS1=MS1/2
      MS2=MS2/2
C..
      IF(MS1.NE.MS2) THEN
        IF(IFLG.EQ.1) call stop_all(t_r, ' PROBLEM IN TWO DIFFER ')
        IF(IFLG.EQ.2) call stop_all(t_r, ' PROBLEM IN ONE DIFFER ')
        IF(IFLG.EQ.3) call stop_all(t_r, ' PROBLEM IN SAME ')
      ENDIF 
C..
      RETURN
      END
