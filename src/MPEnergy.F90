! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
!See 5/1/07
SUBROUTINE AddMPEnergy(Hij,iV,iMaxOrder,Arr,nBasis,iPath,nEl,tLog,ECore,MPEs)
   use constants, only: dp
   use CalcData , only : TMODMPTHEORY,TENPT
   IMPLICIT NONE
   INTEGER iV,nEl,nBasis,iMaxOrder
   HElement_t Hij(0:iV,0:iV)
   HElement_t V(1:iV,1:iV)
   real(dp) Arr(nBasis,2)
   INTEGER iPath(nEl,0:iV)
   LOGICAL tLog,tLogged
   INTEGER i,j
   HElement_t Fi(1:iV),E1
   INTEGER iOrder
   real(dp) MPEs(2:iMaxOrder),E,ECore
   HElement_t MPE
!E1 is the HF Energy.  Ei are Fock energy differences.
   MPE=ECore
   DO i=1,iV
!      EX(1,1)=nEl
!      CALL GetExcitation(iPath(1,0),iPath(i,0),nEl,EX,tSign)
      IF(tENPT) THEN
!Use Epstein-Nesbet denominator
         Fi(i)=Hij(i-1,i-1)
      ELSE
!Use standard MP denominator
         Fi(i)=ECore
         DO j=1,nEl
            Fi(i)=Fi(i)+(Arr(iPath(j,i-1),2))
         ENDDO
         IF(tModMPTheory.and.i.gt.1) THEN
!Encoded in the diagonal part of Hij is the <ij||ij>+<ab||ab> sum for the modified MP Theory
!            call writepathex(6,iPath,iV,nEl,.true.)
!            write(6,*) fi(i)-fi(1),fi(i)+hij(i-1,i-1)-fi(1), hij(i-1,i-1)
            Fi(i)=fi(i)+Hij(i-1,i-1)
         ENDIF
      ENDIF
   ENDDO
   tLogged=.FALSE.
   E1=Hij(0,0)
   E1=E1-Fi(1)
   CALL CalcVij(Hij,V,Fi,E1,iV)
   DO iOrder=2,iMaxOrder
      MPE=0.0_dp
      SELECT CASE(iV)
      CASE(2)
         IF(iOrder.EQ.2) MPE=MPE+V(1,2)*V(2,1)/(Fi(1)-Fi(2))
         IF(iOrder.EQ.3) MPE=MPE+V(1,2)*V(2,2)*V(2,1)/((Fi(1)-Fi(2))**2.0_dp)
         IF(iOrder.EQ.4) MPE=MPE+V(1,2)*V(2,1)*(-V(2,1)*V(1,2)+V(2,2)**2.0_dp)/((Fi(1)-Fi(2))**3.0_dp)
      CASE(3)
         IF(iOrder.EQ.3) MPE=MPE+(V(1,2)*V(2,3)*V(3,1)+V(1,3)*V(3,2)*V(2,1)) &
                                             /((Fi(1)-Fi(2))*(Fi(1)-Fi(3)))
         IF(iOrder.EQ.4) MPE=MPE+ &
                            ((V(1,2)*V(2,3)*V(3,2)*V(2,1)+ &
                             V(1,2)*V(2,2)*V(2,3)*V(3,1)+ &
                             V(1,3)*V(3,2)*V(2,2)*V(2,1)-  &
                             V(1,2)*V(2,1)*V(1,3)*V(3,1))/(Fi(1)-Fi(2)) + &
                            (V(1,2)*V(2,3)*V(3,3)*V(3,1)+ &
                             V(1,3)*V(3,3)*V(3,2)*V(2,1)+ &
                             V(1,3)*V(3,2)*V(2,3)*V(3,1)-  &
                             V(1,2)*V(2,1)*V(1,3)*V(3,1))/(Fi(1)-Fi(3))) &
                           /((Fi(1)-Fi(2))*(Fi(1)-Fi(3)))
      CASE(4)
         IF(iOrder.EQ.4) MPE=MPE+(V(1,2)*V(2,3)*V(3,4)*V(4,1)+ &
                             V(1,2)*V(2,4)*V(4,3)*V(3,1)+ &
                             V(1,3)*V(3,2)*V(2,4)*V(4,1)+ &
                             V(1,3)*V(3,4)*V(4,2)*V(2,1)+ &
                             V(1,4)*V(4,2)*V(2,3)*V(3,1)+ &
                             V(1,4)*V(4,3)*V(3,2)*V(2,1))  &
                           /((Fi(1)-Fi(2))*(Fi(1)-Fi(3))*(Fi(1)-Fi(4)))

      END SELECT
      E=MPE
      MPEs(iOrder)=MPEs(iOrder)+E
      IF(TLOG.AND. ABS(E) .GT. 1.0e-9_dp) THEN
         IF(iOrder.EQ.2) CALL WRITEPATH(13,IPATH,2,NEL,.FALSE.)
         WRITE(13,"(G25.16)",advance='no') E
         tLogged=.TRUE.
      ENDIF
   ENDDO
   IF(tLogged) WRITE(13,*)
   RETURN
END

! V=H1-E1
! H1=H-F
! E1=<Di|H1|Di>
! Vij=<Di|V|Dj>
! Fj=<Dj|F|Dj>
SUBROUTINE CalcVij(Hij,Vij,Fi,E1,iV)
   use constants, only: dp
   IMPLICIT NONE
   INTEGER iV,i,j
   HElement_t Hij(1:iV+1,1:iV+1)
   HElement_t Vij(1:iV,1:iV)
   HElement_t Fi(1:iV),E1
   DO i=1,iV
      DO j=1,iV
         Vij(i,j)=Hij(i,j)
         IF(i.EQ.j) Vij(i,j)=Vij(i,j)-Fi(i)-E1
      ENDDO
   ENDDO
END


!.. Calculate the contribution to the MP2 energy from 
!.. the determinant making this a 2-v graph.
      SUBROUTINE ADDMP2E(HIJS,ARR,NBASIS,IPATH,NEL,TLOG,MP2E)
         use constants, only: dp
         use CalcData , only : TLADDER
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NBASIS
         HElement_t HIJS(0:2)
         real(dp) ARR(NBASIS,2)
         INTEGER IPATH(NEL,0:2)
         INTEGER NI(NEL),NJ(NEL)
         real(dp) MP2E
         real(dp) DENOM,CONTR
         INTEGER I,J,S
         LOGICAL TLOG
         LOGICAL iscsf_old


!.. If we have CSFs, unCSF the elecs
         IF(iscsf_old(IPATH(1,0),NEL)) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(IPATH(I,0),NI(I),S)
            ENDDO
         ELSE
            CALL NECI_ICOPY(NEL,IPATH(1,0),1,NI,1)
         ENDIF
         IF(iscsf_old(IPATH(1,1),NEL)) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(IPATH(I,1),NJ(I),S)
            ENDDO
         ELSE
            CALL NECI_ICOPY(NEL,IPATH(1,1),1,NJ,1)
         ENDIF
!         CALL WRITEDET(13,NI,NEL,.FALSE.)
!         CALL WRITEDET(13,NJ,NEL,.TRUE.)
!.. First find which orbitals have been excited
         I=1
         J=1
         DENOM=0.0_dp
         DO WHILE (I.LE.NEL.OR.J.LE.NEL)
!            WRITE(13,*) I,J,NI(I),NJ(J)
            IF(J.GT.NEL) THEN
!.. 0 has a lower orb than 1, so 0's must have been removed to get 1
!               WRITE(13,*) NI(I),-ARR(NI(I),2)
               DENOM=DENOM-ARR(NI(I),2)
               I=I+1
            ELSEIF(I.GT.NEL) THEN
!.. 1 has a lower orb than 0, so 1's must have been added to get 1
!               WRITE(13,*) NJ(J),ARR(NJ(J),2)
               DENOM=DENOM+ARR(NJ(J),2)
               J=J+1
            ELSEIF(I.LE.NEL.AND.NI(I).LT.NJ(J)) THEN
!.. 0 has a lower orb than 1, so 0's must have been removed to get 1
!               WRITE(13,*) NI(I),-ARR(NI(I),2)
               DENOM=DENOM-ARR(NI(I),2)
               I=I+1
            ELSEIF(J.LE.NEL.AND.NI(I).GT.NJ(J)) THEN
!.. 1 has a lower orb than 0, so 1's must have been added to get 1
!               WRITE(13,*) NJ(J),ARR(NJ(J),2)
               DENOM=DENOM+ARR(NJ(J),2)
               J=J+1
            ELSE
               I=I+1
               J=J+1
            ENDIF
         ENDDO
         IF(tLadder) then
!Ladder sum perturbation theory
            DENOM=DENOM+abs(HIJS(1))**2
         ENDIF
         CONTR=abs(HIJS(1))**2/DENOM
         IF(TLOG.AND.CONTR.GT.1.0e-9_dp) THEN
            CALL WRITEPATH(13,IPATH,2,NEL,.FALSE.)
            WRITE(13,"(G25.16)",advance='no') -CONTR
            WRITE(13,*) HIJS(1),DENOM
         ENDIF
         MP2E=MP2E-CONTR
         RETURN
      END

      Subroutine ModMPDiagElement(hEl,nI,nJ,nEl)
         use Integrals_neci, only : GetUMatEl
         use constants, only: dp
         use SystemData, only: BasisFN
         implicit none
         HElement_t hEl
         integer nEl,nI(nEl),nJ(nEl)
         integer Ex(2,2),ex2(2,2)
         logical tSign
         Ex(1,1)=2
         Call GetExcitation(nI,nJ,nEl,EX,tSign)
         ex2=ex
         ex=ex+1
         ex=ex/2
         if(Ex(1,2).gt.0) then
            hEl= GetUMatEl(ex(1,1),ex(1,2),ex(1,1),ex(1,2))
            if(mod(ex2(1,1)+ex2(1,2),2).eq.0)                                                         & 
     &         hEl=hEl-GetUMatEl(ex(1,1),ex(1,2),ex(1,2),ex(1,1))
            hEl=hEl+GetUMatEl(ex(2,1),ex(2,2),ex(2,1),ex(2,2))
            if(mod(ex2(2,1)+ex2(2,2),2).eq.0)                                                         &
     &         hEl=hEl-GetUMatEl(ex(2,1),ex(2,2),ex(2,2),ex(2,1))
         else
            hEl=0.0_dp
         endif
      End Subroutine
