!Diagonalize a compressed hamiltonian (in HAMIL), returning eigenvalues in W and eigenvectors in CK.
       SUBROUTINE HDIAG(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,LENHAMIL,
     &                     NBLOCKSTARTS,NBLOCKS,BLOCKSYM)
         USE HElem
         IMPLICIT NONE
         INTEGER NDET,LENHAMIL
         INTEGER NROW(NDET),LAB(*)
         TYPE(HElement)  HAMIL(*),CK(NDET,NDET),WORK2,WORK
         REAL*8 W(NDET)
         INTEGER I,J,IND,INDZ,INFO,LDA,NBS
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS),ISUB
         INCLUDE 'basis.inc'
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
C.. CK now contains the eigenvectors, and W the eigenvalues
         WRITE(6,"(A,F19.11,I4)") "GROUND E=",GSEN
         CALL TIHALT('HDIAG     ',ISUB)
         RETURN
      END
 
      SUBROUTINE HDIAG_NH(NDET,NBLOCKSTARTS,NBLOCKS,
     &            NEL,NMRKS,NBASISMAX,NBASIS,G1,NMSH,BRR,
     &            FCK,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,
     &            BETA,I_P,ILOGGING,IFDET,ARR,BLOCKSYM)
         IMPLICIT NONE
         INTEGER NDET,NEL,NBLOCKS,NBLOCKSTARTS(NBLOCKS),NMRKS(NEL,NDET)
         INTEGER NBASISMAX(5,2),NBASIS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) ,G1(*)
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
               CALL WRITEDET(15,NMRKS(1,NBS+III-1),NEL,.FALSE.)
               WRITE(15,"(I7,F19.11,I5)") NBS+III-1,W(III),IDEG

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
