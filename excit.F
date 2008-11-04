!This is an O[2*N] operation for doubles or O[N] operation for single excitations.
!It takes a determinant nI, and a 2x2 matrix indicating the excited 'to' and 'from'
!orbitals, and finds the determinant nJ, while keeping the final determinant ordered.
!This will save having to order the final determinant.

!ExcitMat(1,*) are the indices in the determinant to vacate from nI (the i,j pair)
!ExcitMat(2,*) are the orbitals to occupy in nJ (the a,b pair) (not the index, but the actual orbital):q
!IC should be 1 or 2, depending on whether it is a double or single excitation
!Single excitations should just have ExcitMat(1,1) and ExcitMat(2,1) with orbital
!information.
!The algorithm could be improved for double excitations by only searching through the
!determinant once, reducing it from an O[2N] to O[N] operation, though would be a little
!more fiddly...
      SUBROUTINE FindExcitDet(ExcitMat,nI,IC)
         use SystemData , only : NEl
         IMPLICIT NONE
         INTEGER :: ExcitMat(2,2),i,j,k
         INTEGER :: nI(NEl),IC
         
!         WRITE(19,*) ExcitMat(:,:),nI
!         IF((IC.gt.2).or.(IC.lt.0)) THEN
!             CALL Stop_All("FindExcitDet","Trying to find an excitation"&
!     &          //" which is more than a double away.")
!         ELSEIF(IC.eq.0) THEN
!Excitation is the same determinant - return
!             RETURN
!         ENDIF

         do k=1,IC
!This needs to be done twice for double excitations of nI
             
!             i=1    !Look through the determinant for the orbital to vacate
!             do while(nI(i).ne.ExcitMat(1,k))
!                 i=i+1
!             enddo
             i=ExcitMat(1,k)

!i is now the index of the excited orbital
             IF(ExcitMat(2,k).gt.nI(ExcitMat(1,k))) THEN
!Excited orbital is a higher-indexed orbital than the original

                 IF(i.eq.NEl) THEN
!We want to replace the final orbital with one larger than it - there is no more rearranging to do.
                     nI(i)=ExcitMat(2,k)
                 ELSE
                     j=i+1  !Start looking for insertion point at the index one above the now empty index
                     do while(nI(j).lt.ExcitMat(2,k))
                         nI(j-1)=nI(j)      !Substitute the orbital into the one below it to maintain ordering
                         IF(j.eq.NEl) THEN
!We have reached the end of the determinant list. The new orbital must go here.
                             j=j+1
                             EXIT
                         ENDIF
                         j=j+1
                     enddo
                     nI(j-1)=ExcitMat(2,k)  !Finally found the correct place to insert the new orbital

!Since we have shifted all the orbitals in a certain range in the first insertion, we need to compensate for
!this by changing the insertion point of the next orbital
                     IF((ExcitMat(1,2).gt.ExcitMat(1,1)).and.           &
     &                    (ExcitMat(1,2).lt.j)) THEN
                         ExcitMat(1,2)=ExcitMat(1,2)-1
                     ENDIF

                 ENDIF
             

!             ELSEIF(ExcitMat(2,k).eq.nI(ExcitMat(1,k))) THEN
!We are replacing an orbital with the same orbital. We want to do nothing to the determinant
!                 CYCLE
             ELSE
!orbital to occupy is a lower-indexed orbital than the original
                 IF(i.eq.1) THEN
!We want to replace the final orbital with one smaller than it - there is no more rearranging to do.
                     nI(1)=ExcitMat(2,k)
                 ELSE
                     j=i-1  !Start looking at the index one below the now empty orbital
                     do while(nI(j).gt.ExcitMat(2,k))
                         nI(j+1)=nI(j)
                         IF(j.eq.1) THEN
!We have reached the beginning of the determinant list. The new orbital must go here.
                             j=j-1
                             EXIT
                         ENDIF
                         j=j-1
                     enddo
                     nI(j+1)=ExcitMat(2,k)  !Found the correct place for the newly-occupied orbital to keep order

!Since we have shifted all the orbitals in a certain range in the first insertion, we need to compensate for
!this by changing the insertion point of the next orbital
                     IF((ExcitMat(1,2).lt.ExcitMat(1,1)).and.           &
     &                   (ExcitMat(1,2).gt.j)) THEN
                         ExcitMat(1,2)=ExcitMat(1,2)+1
                     ENDIF

                 ENDIF

             ENDIF

         enddo  !End loop over number of excitations

!         WRITE(19,*) nI(:)

         RETURN
     
      END SUBROUTINE FindExcitDet


      LOGICAL FUNCTION AREDETSEXCITS(D1,D2,NEL,NBASISMAX,ICE)
         IMPLICIT NONE
         INTEGER ICE,NEL,D1(NEL),D2(NEL),NBASISMAX(*)
         LOGICAL ISUHFDET
         INTEGER IGETEXCITLEVEL_
         ICE=IGETEXCITLEVEL_(D1,D2,NEL)
         IF(ICE.LE.2) THEN
C.. The new node is connected to this node
            IF(ISUHFDET(D1,NEL).OR.
     &         ISUHFDET(D2,NEL)) THEN
               IF(ICE.NE.1) THEN
                  AREDETSEXCITS=.TRUE.
               ELSE
                  AREDETSEXCITS=.FALSE.
               ENDIF
            ELSE
               AREDETSEXCITS=.TRUE.
            ENDIF
         ELSE
            AREDETSEXCITS=.FALSE.
         ENDIF
         RETURN
      END

      LOGICAL FUNCTION ISVALIDDET(DET,NEL)
         IMPLICIT NONE
         INTEGER NEL,DET(NEL)
         INTEGER I
         IF(DET(1).LT.1) THEN
            ISVALIDDET=.FALSE.
            RETURN
         ENDIF
         DO I=2,NEL
            IF(DET(I-1).GE.DET(I)) THEN
               ISVALIDDET=.FALSE.
               RETURN
            ENDIF
         ENDDO
         ISVALIDDET=.TRUE.
         RETURN
      END

C.. Get the orbitals which are excited in going from I to J
C.. EX(1,*) are in I, and EX(2,*) are in J
C.. TSIGN is set to the parity of the permutations required to line up the orbitals
C..   TRUE means ODD.
C.. If there are too many excitations to fit, then we put -excitlevel in EX(1,1) and EX(2,1)
!  EX(1,1) is the max number of excitations (passed in as a parameter)
      SUBROUTINE GETEXCITATION(NI,NJ,NEL,EX,TSIGN)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),EX(2,*)
         INTEGER I,J,IPAR
         INTEGER IC1,IC2
         INTEGER iMaxExcit
         LOGICAL TSIGN
         iMaxExcit=EX(1,1)
         EX(1:2,1:iMaxExcit)=0
         IC1=0
         IC2=0
         I=1
         J=1
         IPAR=0
C         CALL WRITEDET(6,NI,NEL,.TRUE.)
C         CALL WRITEDET(6,NJ,NEL,.TRUE.)
         DO WHILE(I.LE.NEL.AND.J.LE.NEL)
C.. Differences from I to J
C            WRITE(6,*) "GE",I,J
            DO WHILE(I.LE.NEL.AND.NI(I).LT.NJ(J))
               IC1=IC1+1
               IF(IC1.LE.iMaxExcit) THEN
                  EX(1,IC1)=NI(I)
                  IPAR=IPAR+I
               ENDIF
               I=I+1
            ENDDO
C.. Differences from J to I
            DO WHILE(I.LE.NEL.AND.J.LE.NEL.AND.NI(I).GT.NJ(J))
               IC2=IC2+1
               IF(IC2.LE.iMaxExcit) THEN
                  EX(2,IC2)=NJ(J)
                  IPAR=IPAR+J
               ENDIF
               J=J+1
            ENDDO
            IF(I.LE.NEL.AND.J.LE.NEL.AND.NI(I).EQ.NJ(J)) THEN
               I=I+1
               J=J+1
            ENDIF
         ENDDO
C.. Deal with remaining I
         DO WHILE(I.LE.NEL)
            IC1=IC1+1
            IF(IC1.LE.iMaxExcit) THEN
               IPAR=IPAR+I
               EX(1,IC1)=NI(I)
            ENDIF
            I=I+1
         ENDDO
C.. Deal with remaining J
         DO WHILE(J.LE.NEL)
            IC2=IC2+1
            IF(IC2.LE.iMaxExcit) THEN
               EX(2,IC2)=NJ(J)
               IPAR=IPAR+J
            ENDIF
            J=J+1
         ENDDO
         IF(iC1.GT.iMaxExcit) THEN
C.. we actually needed more space.  Just list the excitation counts (-ve)
            DO i=1,iMaxExcit
               IF (i.EQ.1) THEN
                  EX(1,1)=-iC1
                  EX(2,1)=-iC2
               ELSE
                  EX(1,i)=0
                  EX(2,i)=0
               ENDIF
            ENDDO
         ELSEIF(iC1.EQ.0) THEN
            EX(1,1)=0
            EX(2,1)=0
         ENDIF
         TSIGN=BTEST(IPAR,0)
         RETURN
      END


C.. Generate excitations of order at most NORDER from NI (excluding NI
Citself) in order
C.. ICLIST is the list of orders of the excitations
C.. NMIN is the minimum level of excitation
C.. IF NORDER=0, we SET NORDER=NEL

C.. THis does not use the excitation generators, but enumerates all possible
C.. determinants within double excitations, which will end up fairly inefficient.
      SUBROUTINE GENEXCIT(NI,NORDER,NBASIS,NEL,LIST,ICLIST,NLIST,NMIN,
     &            G1,TSYM,NBASISMAX,TCOUNT)
         use CalcData , only : TNEWEXCITATIONS
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NORDER,NBASIS,LIST(NEL,*),ICLIST(*)
         INTEGER NEXCIT(NEL),NLIST,nBasisMax(5,*)
         TYPE(BASISFN) ISYM
         INTEGER NO,IFDET
         LOGICAL TSYM,TCOUNT
         INTEGER G1(*),NMIN

         INTEGER*8 STORE(6)
         INTEGER ICOUNT,ILEVEL


         NO=NORDER
         IF(NO.EQ.0) NO=NEL
         IF(TSYM) CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         IF(NBASISMAX(4,7).EQ.1) THEN
            NLIST=1
C.. If we're actually using CSFs
            CALL GNCSFS(NEL,NBASIS,0,NBASISMAX,
     &         LIST,TCOUNT,G1,TSYM,ISYM%MS,TSYM,ISYM,NLIST,IFDET,
     &         .TRUE.,NI,ICLIST,.TRUE.,NO)
         ELSE
            IF(TNEWEXCITATIONS.AND.NORDER.LE.2) THEN
               IF(NORDER.EQ.1) ILEVEL=1
               IF(NORDER.EQ.2) THEN
                  IF(NMIN.EQ.1) THEN
                     ILEVEL=3
                  ELSE
                     ILEVEL=2
                  ENDIF
               ENDIF
C.. use a new excitation routine - THIS DOES NOT PRODUCE AN ORDERED LIST
               CALL SYMSETUPEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &            TCOUNT,ICOUNT,ILEVEL,1,nEl)
               IF(.NOT.TCOUNT) THEN
                  IF(ICOUNT.GT.NLIST) THEN
                     WRITE(6,*) "Exctiation list size",NLIST
                     WRITE(6,*) "Requires size",ICOUNT
                     STOP "EXCITATIONS LIST NOT LARGE ENOUGH"
                  ENDIF
                  CALL SYMGENEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &               LIST,ICLIST,ICOUNT)
               ENDIF
               NLIST=ICOUNT
            ELSE
               NLIST=1
               CALL GENEXCIT_R(NI,NO,NEL,1,NBASIS,LIST,ICLIST,
     &            NLIST,NEXCIT,NO,1,NMIN,G1,ISYM,TSYM,NBASISMAX,TCOUNT)
               NLIST=NLIST-1
            ENDIF
         ENDIF
         RETURN
      END

C.. AT 29/1/04 - Looks like this won't scale very well to lots of
Celectrons
C..      as there will be very few which are within NORDER excitations
Cof
C..      our original det. (Or for that matter a large basis)
C.. Recursively go through each electron ( we're on NELEC), filling
CLIST.
C.. with up to NORDER excitations.
C.. NEXCIT is the determinant being constructed.
C.. NLEFT is the number of possible excitations left
C.. NSTARTFN is the basis fn to start with for this electron
      RECURSIVE SUBROUTINE GENEXCIT_R(NI,NORDER,NEL,NELEC,NBASIS,LIST,
     &         ICLIST,NLISTPOS,NEXCIT,NLEFT,NSTARTFN,NMIN,G1,ISYM,TSYM,
     &         NBASISMAX,TCOUNT)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NELEC,NORDER,NLISTPOS,NLEFT,NBASIS,NSTARTFN,NNLEFT
         INTEGER NI(NEL),NEXCIT(NEL),LIST(NEL,*),ICLIST(*)
         INTEGER I,J,I1,I2,NMIN,NMAXEX,K
         LOGICAL LISINOLD,LSYM,TSYM,LCHKSYM,TCOUNT
         TYPE(BASISFN) G1(NBASIS),ISYM,ISYM2
         INTEGER nBasisMax(5,*)
C.. NMAXEX is the maximum number of excitations we're allowed to  have
Cleft
         NMAXEX=NORDER-NMIN
C..  I is the new basis fn for electron NELEC
C.. We see if it's in the original determinant
C.. we need to leave at least NEL-NELEC basis fns for the following
Celectrons
         DO I=NSTARTFN,NBASIS-(NEL-NELEC)
            LISINOLD=.FALSE.
            DO J=1,NEL
               IF(NI(J).EQ.I) LISINOLD=.TRUE.
            ENDDO
            NEXCIT(NELEC)=I
            NNLEFT=NLEFT
            IF(.NOT.LISINOLD) NNLEFT=NNLEFT-1
C.. NNLEFT is the number of excitations left
C.. check to see if we're actually exciting an electron, and if we've
Cany
C.. spare excitations left
C.. if we've allocated the last electron, we need to see if the det
C.. is allowed, and if so, store it
            IF(NELEC.EQ.NEL) THEN
C.. if we've excited at least one elec
               IF(NNLEFT.LE.NMAXEX.AND.NNLEFT.GE.0) THEN
C.. we check whether we're allowed this excitation owing to sym
C                  CALL WRITEDET(6,NEXCIT,NEL,.TRUE.)
                  IF(TSYM) THEN
                     CALL GETSYM(NEXCIT,NEL,G1,NBASISMAX,ISYM2)
                     LSYM=LCHKSYM(ISYM,ISYM2)
                  ELSE
                     LSYM=.TRUE.
                  ENDIF
C                  WRITE(6,*) ISYM,ISYM2,LSYM
                  IF(LSYM) THEN
                     IF(.NOT.TCOUNT) THEN
                        CALL NECI_ICOPY(NEL,NEXCIT,1,
     &                                           LIST(1:NEL,NLISTPOS),1)
                        ICLIST(NLISTPOS)=NORDER-NNLEFT
                     ENDIF
                     NLISTPOS=NLISTPOS+1
                  ENDIF
               ENDIF
            ELSEIF(NNLEFT.GE.0) THEN
C.. we start the next electron on the basis fn after this one
               CALL GENEXCIT_R(NI,NORDER,NEL,NELEC+1,NBASIS,LIST,ICLIST,
     &            NLISTPOS,NEXCIT,NNLEFT,I+1,NMIN,G1,ISYM,TSYM,
     &            NBASISMAX,TCOUNT)
            ENDIF
         ENDDO
         RETURN
      END


