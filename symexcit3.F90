MODULE SymExcit3
! This module contains excitation generators able to enumerate all possible excitations given a starting determinant.
! Unlike symexcit.F90 however, these excitation generators are able to deal with cases where the alpha and beta orbitals 
! have different symmetries.  This is particularly relevant when dealing with certain unrestricted cases, or when we
! are truncating (or freezing) orbitals in such a way as to remove different alpha symm irreps from the beta.

    USE SystemData, only: NEl,G1,NIfD,nBasis,tNoSymGenRandExcits
    USE GenRandSymExcitNUMod, only: SymLabelList2,SymLabelCounts2,ClassCountInd,ScratchSize
    IMPLICIT NONE


    CONTAINS


    SUBROUTINE CountExcitations3(nI,exflag,nSingleExcits,nDoubleExcits)
! This routine simply counts the excitations in terms of single and doubles from the nI determinant.    
! The exflag sent through indicates which should be counted - exflag=1 means only singles, exflag=2 means
! only doubles, and anything else both are counted.
        USE SymData, only: nSymLabels
        USE SystemData , only: ElecPairs
        USE GenRandSymExcitNUMod , only: PickElecPair,ConstructClassCounts,ClassCountInd,ScratchSize 
        INTEGER :: nSingleExcits,nDoubleExcits,Symi,i,j,Spini,nI(NEl)
        INTEGER :: iSpn,Elec1Ind,Elec2Ind,SymProduct,exflag
        INTEGER :: Syma,Symb,Spina,Spinb,StartSpin,EndSpin
        INTEGER :: ClassCount2(ScratchSize),SumMl
        INTEGER :: ClassCountUnocc2(ScratchSize)

        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
! This sets up arrays containing the number of occupied and unoccupied in each symmetry.
! ClassCounts2(1,:)=No alpha occupied, ClassCounts2(2,:)=No Beta occupied.
! ClassCountsUnocc2(1,:)=No alpha unocc, ClassCounts2Unocc2(2,:)=No Beta unocc.
! The second index of these arrays referrs to the symmetry (0 -> 7).

! Only counting.  Run through each occupied orbital, and count the number of spin and symmetry allowed orbitals it 
! may be excited to.
        nSingleExcits=0
        nDoubleExcits=0

        IF(exflag.ne.2) THEN
! Count the singles.            
! Take each electron and find out the number of symmetry allowed orbitals it may be excited to.
            do i=1,NEl
                IF(tNoSymGenRandExcits) THEN
                    Symi=0
                ELSE
                    Symi=INT(G1(nI(i))%Sym%S,4)
                ENDIF
                IF((G1(nI(i))%Ms).eq.-1) Spini=2        ! G1(i)%Ms is -1 for beta, and 1 for alpha.
                IF((G1(nI(i))%Ms).eq.1) Spini=1         ! Translate this into 1 for alpha and 2 for beta
                                                        ! for the ClassCount arrays.

! This electron in orbital of SymI and SpinI can only be excited to orbitals with the same spin and symmetry.                
! Then add in the number of unoccupied orbitals with the same spin and symmetry to which each electron may be excited.
            
                nSingleExcits=nSingleExcits+ClassCountUnocc2(ClassCountInd(Spini,Symi,0))

            enddo
        ENDIF

! This is the end of the singles.        
!            WRITE(6,*) 'Number of singles',nSingleExcits

! For the doubles, first pick an electron pair i,j.
! Based on these orbitals, run through each spin and each symmetry - take this to be orbital a.
! Multiply the number with these symmetries by the number of possible b orbitals which correspond.
! Do this for all a and then all i,j pairs.
        IF(exflag.ne.1) THEN
            do i=1,ElecPairs

! iSpn=2 for alpha beta pair, ispn=3 for alpha alpha pair and ispn=1 for beta beta pair.
                CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,i)

                StartSpin=1
                EndSpin=2
                IF(iSpn.eq.3) EndSpin=1
                IF(iSpn.eq.1) StartSpin=2
                do Spina=StartSpin,EndSpin            ! Run through both spins, orbital a may be alpha or beta.
                    IF(iSpn.eq.2) THEN
! Spin of orbital b should be opposite to orbital a.                    
                        IF(Spina.eq.1) Spinb=2
                        IF(Spina.eq.2) Spinb=1
                    ELSE
! Spin of orbital b should be opposite to orbital a.                    
                        IF(Spina.eq.1) Spinb=1
                        IF(Spina.eq.2) Spinb=2
                    ENDIF

                    do Syma=0,nSymLabels-1

! Need to work out the symmetry of b, given the symmetry of a (Sym).                    
                        Symb=IEOR(Syma,SymProduct)

                        IF((Spina.eq.Spinb).and.(Syma.eq.Symb)) THEN
                            ! If the spin and spatial symmetries of a and b are the same
                            ! there will exist a case where Orba = Orbb, want to remove this.
                            nDoubleExcits=nDoubleExcits+(ClassCountUnocc2(ClassCountInd(Spina,Syma,0))*(ClassCountUnocc2(ClassCountInd(Spinb,Symb,0))-1))
                        ELSE
                            nDoubleExcits=nDoubleExcits+(ClassCountUnocc2(ClassCountInd(Spina,Syma,0))*ClassCountUnocc2(ClassCountInd(Spinb,Symb,0)))
                        ENDIF
                    enddo

                enddo
            enddo
            nDoubleExcits=nDoubleExcits/2

        ENDIF

!        WRITE(6,*) 'Number of doubles',nDoubleExcits

    ENDSUBROUTINE CountExcitations3



    SUBROUTINE GenExcitations3(nI,iLut,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
! This routine finds in turn, every possible excitation from determinant nI.
! The excited determinant is then returned as nJ.
! exflag indicates which excitations we want to find.  exflag=1 - only singles are returned, exflag=2 - only
! doubles are returned and anything else returns the singles followed by the doubles.
! ExcitMat3 holds the orbitals involved in the excitation.
! If an excitation matrix of 0's is passed through, the first single or double is found.
! After this, the routine reads in the ExcitMat and finds the next excitation after this.
! ExcitMat(1,*) are the **indices** in the determinant to vacate from nI (the i,j pair)
! ExcitMat(2,*) are the orbitals to occupy in nJ (the a,b pair) (not the index, but the actual orbital)
! If tParity is true, two orbitals need to be switched in order to better represent the excitation, therefore a 
! negative sign must be included when finding the H element.
! When there are no more symmetry allowed excitations, tAllExcitFound becomes true.
        INTEGER :: nI(NEl),iLut(0:NIfD),nJ(NEl),nSingles,nDoubles,ExcitMat3(2,2),exflag
        LOGICAL :: tCountOnly,tAllExcitFound,tParity

        IF(exflag.eq.2) THEN
! Just generate doubles            

            CALL GenDoubleExcit(nI,iLut,nJ,ExcitMat3,tParity,tAllExcitFound)

        
        ELSE 
! Generate singles, returning Orbi and Orba as non-zero, but keeping the others 0.        

            CALL GenSingleExcit(nI,iLut,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
            
            ! When the last single is input, providing exflag is not 1, the first double is then found
            ! and from then on GenDoubleExcit is called.

        ENDIF


    ENDSUBROUTINE GenExcitations3



    SUBROUTINE GenSingleExcit(nI,iLut,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
! Despite being fed four indices, this routine finds single excitations.  Orbi -> Orba. (Orbj and Orbb remain 0).
! Feeding in 0 indices indicates it is the first excitation that needs to be found.
! The single excitation goes from orbital i to a, from determinant nI to nJ.
! When the last single is found it then finds the first double excitation, unless exflag=1 in which tAllExcitFound 
! becomes true and no more excitations are generated.
        USE SymData, only: nSymLabels
        INTEGER :: i,a,nI(NEl),Orbi,Orba,Symi,Finala,iLut(0:NIfD),nJ(NEl)
        INTEGER :: Orbj,Orbb,NoOcc,k,ExcitMat3(2,2),exflag,SymInd
        LOGICAL :: tInitOrbsFound,tParity,tAllExcitFound
        INTEGER , SAVE :: OrbiIndex,OrbaIndex,Spini,NewSym,OldSym

!        WRITE(6,*) 'Original Determinant',nI

        tInitOrbsFound=.false.
        Orbi=ExcitMat3(1,1)
        Orba=ExcitMat3(2,1)

        IF((Orbi.eq.0).or.(Orba.eq.0)) THEN           ! Want to find the first excitation.

            OrbiIndex=1
            Orbi=nI(OrbiIndex)                              ! Take the first occupied orbital

            IF(tNoSymGenRandExcits) THEN
                Symi=0
            ELSE
                Symi=INT(G1(Orbi)%Sym%S,4)                      ! and find its spin and spat symmetries.
            ENDIF
            IF((G1(Orbi)%Ms).eq.-1) Spini=2  
            IF((G1(Orbi)%Ms).eq.1) Spini=1  
            OrbaIndex=SymLabelCounts2(1,ClassCountInd(Spini,Symi,0))  ! Start considering a at the first allowed symmetry.

        ELSE
            Orbi=nI(OrbiIndex)                              ! Begin by using the same i as last time - check if there are any 
                                                            ! more possible excitations from this.

! At this stage, OrbaIndex is the a from the previous excitation.
            SymInd=ClassCountInd(Spini,INT(G1(Orbi)%Sym%S,4),0)

            IF(OrbaIndex.eq.(SymLabelCounts2(1,SymInd)+SymLabelCounts2(2,SymInd)-1)) THEN
                !Orba was the last in the symmetry block. Do not allow OrbaIndex+1

! Either we're got to the final spin symmetry, or the next orbital after Orba does not have the same symmetry as Orbi.                
! Need to move onto the next i, and find a new a to match.
                OrbiIndex=OrbiIndex+1
                IF(OrbiIndex.le.NEl) THEN
                    Orbi=nI(OrbiIndex)
                    IF(tNoSymGenRandExcits) THEN
                        Symi=0
                    ELSE
                        Symi=INT(G1(Orbi)%Sym%S,4)                  
                    ENDIF
                    IF((G1(Orbi)%Ms).eq.-1) Spini=2  
                    IF((G1(Orbi)%Ms).eq.1) Spini=1  
                    OrbaIndex=SymLabelCounts2(1,ClassCountInd(Spini,Symi,0))
                ENDIF

            ELSE
! There are more possible excitations from orbital a, simply check the next orbital after the current a.
                OrbaIndex=OrbaIndex+1

                IF(tNoSymGenRandExcits) THEN
                    Symi=0
                ELSE
                    Symi=INT(G1(Orbi)%Sym%S,4)           
                ENDIF
            ENDIF
        ENDIF

        do while (.not.tInitOrbsFound)

            IF(OrbiIndex.gt.NEl) THEN
! If we've read in the last single, set orbi, orbj, orba, and orbb to 0 and call gendoubleexcit.        
                IF(exflag.ne.1) THEN
                    ExcitMat3(:,:)=0
                    CALL GenDoubleExcit(nI,iLut,nJ,ExcitMat3,tParity,tAllExcitFound)
                    exflag=2
                ELSE
                    tAllExcitFound=.true.
                ENDIF
                EXIT
            ENDIF

! To find Orba, take the first in SymLabelList2 with the same symmetry and spin.                
! SymLabelCounts2(spin,1,symmetry) gives the index in SymLabelList2 where that spin and symmetry starts.                
            Orba=SymLabelList2(OrbaIndex)

            SymInd=ClassCountInd(Spini,INT(G1(Orbi)%Sym%S,4),0)

! Need to also make sure orbital a is unoccupied, so make sure the orbital is not in nI.
            NoOcc=0
            do while (BTEST(iLut((Orba-1)/32),MOD((Orba-1),32))) 
! While this is true, Orba is occupied, so keep incrementing Orba until it is not.                    
                NoOcc=NoOcc+1
                Orba=SymLabelList2(OrbaIndex+NoOcc)
                IF((OrbaIndex+NoOcc).gt.(SymLabelCounts2(1,SymInd)+SymLabelCounts2(2,SymInd)-1)) EXIT
            enddo

! Then check we have not overrun the symmetry block while skipping the occupied orbitals.                
            IF(tNoSymGenRandExcits) THEN
                NewSym=0
            ELSE
                NewSym=INT(G1(Orba)%Sym%S,4)
            ENDIF
            IF(NewSym.eq.Symi) THEN
! If not, then these are the new Orbi and Orba.                
                tInitOrbsFound=.true.
                OrbaIndex=OrbaIndex+NoOcc
            ELSE

! If we have, move onto the next occupied orbital i, no symmetry allowed single excitations exist from the first.                
                OrbiIndex=OrbiIndex+1
                IF(OrbiIndex.le.NEl) THEN
                    Orbi=nI(OrbiIndex)
                    IF(tNoSymGenRandExcits) THEN
                        Symi=0
                    ELSE
                        Symi=INT(G1(Orbi)%Sym%S,4)                      ! and find its spin and spat symmetries.
                    ENDIF
                    IF((G1(Orbi)%Ms).eq.-1) Spini=2  
                    IF((G1(Orbi)%Ms).eq.1) Spini=1  
                    OrbaIndex=SymLabelCounts2(1,ClassCountInd(Spini,Symi,0))
                ENDIF
            ENDIF

        enddo

        IF(ExcitMat3(1,2).eq.0) CALL FindNewSingDet(nI,nJ,OrbiIndex,OrbA,ExcitMat3,tParity)
            
!        WRITE(6,*) 'Excitation is from :',ExcitMat3(1,1),' to ',ExcitMat3(2,1)
!        CALL FLUSH(6)
!        WRITE(6,*) 'The new determinant is :',nJ(:)


    ENDSUBROUTINE GenSingleExcit



    SUBROUTINE GenDoubleExcit(nI,iLut,nJ,ExcitMat3,tParity,tAllExcitFound)
! This generates one by one, all possible double excitations.
! This involves a way of ordering the electron pairs i,j and a,b so that given an i,j and a,b we can find the next.
! The overall symmetry must also be maintained - i.e. if i and j are alpha and beta, a and b must be alpha and beta
! or vice versa.
        USE SystemData , only: ElecPairs
        USE GenRandSymExcitNUMod , only: PickElecPair,FindNewDet 
        INTEGER :: nI(NEl),iLut(0:NIfD),Orbj,Orbi,Orba,Orbb,OrbbSpin,Syma,Symb,NewSym,SymInd
        INTEGER :: Elec1Ind,Elec2Ind,SymProduct,iSpn,Spinb,nJ(NEl),i,k,ExcitMat3(2,2),SumMl
        INTEGER , SAVE :: ijInd,OrbaIndex,OrbbIndex,Spina
        LOGICAL :: tDoubleExcitFound,tFirsta,tFirstb,tNewij,tNewa,tAllExcitFound,tParity

        tDoubleExcitFound=.false.
        tFirsta=.false.
        tFirstb=.false.

        Orbi=ExcitMat3(1,1)
        Orbj=ExcitMat3(1,2)
        Orba=ExcitMat3(2,1)
        Orbb=ExcitMat3(2,2)
!        WRITE(6,*) 'Orbi,Orbj,Orba,Orbb',Orbi,Orbj,Orba,Orbb
!        CALL FLUSH(6)

        do while (.not.tDoubleExcitFound)

            IF(Orbi.eq.0) THEN
                ijInd=1
! If Orbi, then we are choosing the first double.             
! It is therefore also the first set of a and b for this electron pair i,j.
                tFirsta=.true.
                tFirstb=.true.
            ENDIF

! Otherwise we use the previous ijInd and the saved indexes for a and b.
! This routine allows us to pick an electron pair i,j specified by the index ijInd.
! The i and j orbitals are then given by nI(Elec1Ind) and nI(Elec2Ind), and the symmetry product of the two is 
! SymProduct and the spin iSpn.
! iSpn=2 for alpha beta pair, ispn=3 for alpha alpha pair and ispn=1 for beta beta pair.
            CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,ijInd)

            tNewij=.false.
! This becomes true when we can no longer find an allowed orbital a for this ij pair and we need to move onto the next.
            do while ((.not.tNewij).and.(.not.tDoubleExcitFound))                  ! This loop runs through the allowed a orbitals
                                                                                   ! until a double excitation is found.


                IF(tFirsta) THEN                    
! If this is the first double we are picking with this ij, we start with the alpha spin, unless i and j are both beta.
! There is no restriction on the symmetry for orbital a - although clearly the symmetry we pick determins b.
                    IF(iSpn.eq.1) THEN
                        Spina=2
                        OrbaIndex=2
                    ELSE
                        Spina=1
                        OrbaIndex=1
                    ENDIF
                ENDIF

! If it is not the first, we have stored the previous spina and orba index - need to start with these and see                
! if any more double remain.
                Orba=SymLabelList2(OrbaIndex)

! The orbital chosen must be unoccupied.  This is just a test to make sure this is the case.
                do while (BTEST(iLut((Orba-1)/32),MOD((Orba-1),32))) 

! If not, we move onto the next orbital.                    
                    IF(iSpn.ne.2) THEN
!Increment by two, since we want to look at the same spin state.
                        OrbaIndex=OrbaIndex+2
                    ELSE
!Increment by one, since we want to look at both alpha and beta spins.
                        OrbaIndex=OrbaIndex+1
                    ENDIF

                    IF(OrbaIndex.gt.ScratchSize) THEN
!We have reached the end of all allowed symmetries for the a orbital, only taking into account spin symmetry. Choose new ij pair now.
                        tNewij=.true.
                        EXIT
                    ENDIF

! Otherwise the new orbital a is the first unoccupied orbital of allowed symmetry etc.                    
                    Orba=SymLabelList2(OrbaIndex)
                enddo

! If we have got to the end of the a orbitals, and need a new i,j pair, we increment ijInd and check
! this hasn't gone beyond the limits - bail out if we have.
                IF(tNewij) THEN
                    ijInd=ijInd+1
                    IF(ijInd.gt.ElecPairs) THEN
                        tDoubleExcitFound=.true.
                        tAllExcitFound=.true.
                        ! AllExcitFound true indicates there are no more symmetry allowed double excitations.
                    ENDIF
                    EXIT
                ENDIF


                tNewa=.false.
                !Find a b
                do while ((.not.tNewa).and.(.not.tDoubleExcitFound))

! We now have i,j,a and we just need to pick b.
! First find the spin of b.
                    IF(iSpn.eq.1) THEN
                        Spinb=2
                    ELSEIF(iSpn.eq.3) THEN
                        Spinb=1
                    ELSE
                        IF(Spina.eq.1) Spinb=2
                        IF(Spina.eq.2) Spinb=1
                    ENDIF
! Then find the symmetry of b.
                    IF(tNoSymGenRandExcits) THEN
                        Syma=0
                    ELSE
                        Syma=INT(G1(Orba)%Sym%S,4)
                    ENDIF
                    Symb=IEOR(Syma,SymProduct)

! If this is the first time we've picked an orbital b for these i,j and a, begin at the start of the symmetry block.
! Otherwise pick up where we left off last time.
                    IF(tFirstb) THEN
                        OrbbIndex=SymLabelCounts2(1,ClassCountInd(Spinb,Symb,0))
                    ELSE
!Update new orbital b index - we want to keep the same spin, and since the orbitals alternate alpha/beta, we increment by two.
                        OrbbIndex=OrbbIndex+2
                    ENDIF

! If the new b orbital is still within the limits, check it is unoccupied and move onto the next orbital if it is.                    
                    IF(tNoSymGenRandExcits) THEN
                        NewSym=0
                    ELSE
                        NewSym=INT(G1(SymLabelList2(OrbbIndex))%Sym%S,4)
                    ENDIF
                    SymInd=ClassCountInd(Spinb,NewSym,0)

                    IF(OrbbIndex.eq.(SymLabelCounts2(1,SymInd)+SymLabelCounts2(2,SymInd)-1)) THEN
! If we have already gone beyond the symmetry limits by choosing the next b orbital, pick a new a orbital.                        
                        tNewa=.true.
                    ELSE
                        Orbb=SymLabelList2(OrbbIndex)
! Checking the orbital b is unoccupied and > a.                        
                        do while ((BTEST(iLut((Orbb-1)/32),MOD((Orbb-1),32))).or.(Orbb.le.Orba))
                            !Orbital is occupied - try again
                            IF(OrbbIndex.eq.(SymLabelCounts2(1,SymInd)+SymLabelCounts2(2,SymInd)-1)) THEN
                                !Reached end of symmetry block - need new a
                                tNewa=.true.
                                EXIT
                            ENDIF

                            !Update new orbital b index - we want to keep the same spin, and since the orbitals alternate alpha/beta, we increment by two.
                            OrbbIndex=OrbbIndex+2
                            Orbb=SymLabelList2(OrbbIndex)
                        enddo
                    ENDIF

! If we are moving onto the next a orbital, check we don't also need a new ij pair.                    
                    IF(tNewa) THEN
                        IF(iSpn.ne.2) THEN
!Increment by two, since we want to look at the same spin state.
                            OrbaIndex=OrbaIndex+2
                        ELSE
!Increment by one, since we want to look at both alpha and beta spins.
                            OrbaIndex=OrbaIndex+1
                        ENDIF
                        tFirstb=.true.
                        IF(OrbaIndex.gt.ScratchSize) THEN
!We have reached the end of all allowed symmetries for the a orbital, only taking into account spin symmetry. Choose new ij pair now.
                            tNewij=.true.
                            ijInd=ijInd+1
                            IF(ijInd.gt.ElecPairs) THEN
                                tAllExcitFound=.true.
                            ENDIF
                        ENDIF
                    ELSE
!If we don't need a new a, we have found an excitation ij -> ab that is accepted.                          
                        tDoubleExcitFound=.true.
                    ENDIF

                enddo

            enddo

! This is the loop for new ij pairs - if we are choosing a new ij we are automatically choosing a new a and b also.            
            
            tFirsta=.true.
            tFirstb=.true.

        enddo


        CALL FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,Orba,Orbb,ExcitMat3,tParity)

!        WRITE(6,*) 'From',ExcitMat3(1,:)
!        WRITE(6,*) 'To',ExcitMat3(2,:)

!        WRITE(6,*) 'Excitation from : ',ExcitMat3(1,1),ExcitMat3(1,2),' to ',Orba,Orbb
!        WRITE(6,*) 'These have symmetries : ',INT(G1(ExcitMat3(1,1))%Sym%S,4),INT(G1(ExcitMat3(1,2))%Sym%S,4),' to ',INT(G1(Orba)%Sym%S,4),INT(G1(Orbb)%Sym%S,4)
!        WRITE(6,*) 'The new determinant is : ',nJ(:)
!        CALL FLUSH(6)


    ENDSUBROUTINE GenDoubleExcit

!This routine creates the final determinant for a single excitation.
    SUBROUTINE FindNewSingDet(nI,nJ,Elec1Ind,OrbA,ExcitMat3,tParity)
        INTEGER :: nI(NEl),nJ(NEl),Elec1Ind,OrbA,ExcitMat3(2,2)
        LOGICAL :: tParity

!First construct ExcitMat3
        ExcitMat3(1,1)=Elec1Ind
        ExcitMat3(2,1)=OrbA
        ExcitMat3(1,2)=0
        ExcitMat3(2,2)=0
        nJ(:)=nI(:)
        CALL FindExcitDet(ExcitMat3,nJ,1,tParity)

    END SUBROUTINE FindNewSingDet


END MODULE SymExcit3

