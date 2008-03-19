!See 5/1/07
SUBROUTINE AddMPEnergy(Hij,iV,iMaxOrder,Arr,nBasis,iPath,nEl,tLog,ECore,MPEs)
   USE HElem
   USE Calc , only : TMODMPTHEORY,TENPT
   IMPLICIT NONE
   INTEGER iV,nEl,nBasis,iMaxOrder
   TYPE(HElement) Hij(0:iV,0:iV)
   TYPE(HElement) V(1:iV,1:iV)
   REAL*8 Arr(nBasis,2)
   INTEGER iPath(nEl,0:iV)
   LOGICAL tLog,tSign,tLogged
   INTEGER i,j
   TYPE(HElement) Fi(1:iV),E1
   INTEGER iOrder
   TYPE(HDElement) MPEs(2:iV),E,ECore
   TYPE(HElement) MPE
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
            Fi(i)=Fi(i)+HElement(Arr(iPath(j,i-1),2))
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
      MPE=0.D0
      SELECT CASE(iV)
      CASE(2)
         IF(iOrder.EQ.2) MPE=MPE+V(1,2)*V(2,1)/(Fi(1)-Fi(2))
         IF(iOrder.EQ.3) MPE=MPE+V(1,2)*V(2,2)*V(2,1)/((Fi(1)-Fi(2))**2.D0)
         IF(iOrder.EQ.4) MPE=MPE+V(1,2)*V(2,1)*(-V(2,1)*V(1,2)+V(2,2)**2.D0)/((Fi(1)-Fi(2))**3.D0)
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
      IF(TLOG.AND. ABS(E%v) .GT. 1.D-9) THEN
         IF(iOrder.EQ.2) CALL WRITEPATH(13,IPATH,2,NEL,.FALSE.)
         WRITE(13,"(G25.16,$)") E
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
   USE HElem
   IMPLICIT NONE
   TYPE(HElement) Hij(1:iV+1,1:iV+1)
   TYPE(HElement) Vij(1:iV,1:iV)
   TYPE(HElement) Fi(1:iV),E1
   INTEGER iV,i,j
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
         USE HElem
         USE Calc , only : TLADDER
         IMPLICIT NONE
         TYPE(HElement) HIJS(0:2)
         REAL*8 ARR(NBASIS,2)
         INTEGER IPATH(NEL,0:2),NEL,NBASIS
         INTEGER NI(NEL),NJ(NEL)
         REAL*8 MP2E
         REAL*8 DENOM,CONTR
         INTEGER I,J,S
         LOGICAL TLOG
         LOGICAL ISCSF


!.. If we have CSFs, unCSF the elecs
         IF(ISCSF(IPATH(1,0),NEL)) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(IPATH(I,0),NI(I),S)
            ENDDO
         ELSE
            CALL ICOPY(NEL,IPATH(1,0),1,NI,1)
         ENDIF
         IF(ISCSF(IPATH(1,1),NEL)) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(IPATH(I,1),NJ(I),S)
            ENDDO
         ELSE
            CALL ICOPY(NEL,IPATH(1,1),1,NJ,1)
         ENDIF
!         CALL WRITEDET(13,NI,NEL,.FALSE.)
!         CALL WRITEDET(13,NJ,NEL,.TRUE.)
!.. First find which orbitals have been excited
         I=1
         J=1
         DENOM=0.D0
         DO WHILE (I.LE.NEL.OR.J.LE.NEL)
!            WRITE(13,*) I,J,NI(I),NJ(J)
            IF(J.GT.NEL.OR.(I.LE.NEL.AND.NI(I).LT.NJ(J))) THEN
!.. 0 has a lower orb than 1, so 0's must have been removed to get 1
!               WRITE(13,*) NI(I),-ARR(NI(I),2)
               DENOM=DENOM-ARR(NI(I),2)
               I=I+1
            ELSEIF(I.GT.NEL.OR.(J.LE.NEL.AND.NI(I).GT.NJ(J))) THEN
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
            DENOM=DENOM+SQ(HIJS(1))
         ENDIF
         CONTR=SQ(HIJS(1))/DENOM
         IF(TLOG.AND.CONTR.GT.1.D-9) THEN
            CALL WRITEPATH(13,IPATH,2,NEL,.FALSE.)
            WRITE(13,"(G25.16,$)") -CONTR
            WRITE(13,*) HIJS(1),DENOM
         ENDIF
         MP2E=MP2E-CONTR
         RETURN
      END

      Subroutine ModMPDiagElement(hEl,nI,nJ,nEl,nBasisMax,UMat,ALat,nBasis,iss,G1)
         use UMatCache, only : GetUMatEl
         USE HElem
         use System, only: BasisFN
         implicit none
         Type(HElement) hEl,UMat(*)
         integer nEl,nI(nEl),nJ(nEl),nBasis
         integer nBasisMax(5,5)
         type(BasisFn) G1(*)
         integer Ex(2,2),ex2(2,2),iss
         logical tSign
         real*8 ALat(3)
         Ex(1,1)=2
         Call GetExcitation(nI,nJ,nEl,EX,tSign)
         ex2=ex
         ex=ex+1
         ex=ex/2
         if(Ex(1,2).gt.0) then
            hEl= GetUMatEl(NBASISMAX,UMAT,ALAT,nBasis,ISS,G1,ex(1,1),ex(1,2),ex(1,1),ex(1,2))
            if(mod(ex2(1,1)+ex2(1,2),2).eq.0)                                                         & 
     &         hEl=hEl-GetUMatEl(NBASISMAX,UMAT,ALAT,nBasis,ISS,G1,ex(1,1),ex(1,2),ex(1,2),ex(1,1))
            hEl=hEl+GetUMatEl(NBASISMAX,UMAT,ALAT,nBasis,ISS,G1,ex(2,1),ex(2,2),ex(2,1),ex(2,2))
            if(mod(ex2(2,1)+ex2(2,2),2).eq.0)                                                         &
     &         hEl=hEl-GetUMatEl(NBASISMAX,UMAT,ALAT,nBasis,ISS,G1,ex(2,1),ex(2,2),ex(2,2),ex(2,1))
         else
            hEl=0.d0
         endif
      End Subroutine
