MODULE SymExcit3
! This module contains excitation generators able to enumerate all possible excitations given a starting determinant.
! Unlike symexcit.F90, these excitation generators are able to deal with cases where the alpha and beta orbitals 
! have different symmetries.  This is particularly relevant when dealing with certain unrestricted cases, or when we
! are truncating (or freezing) orbitals in such a way as to remove different alpha symm irreps from the beta.


    USE SystemData, only: NEl,G1,NIfD,nBasis 
    USE GenRandSymExcitNUMod, only: SymLabelList2,SymLabelCounts2
    IMPLICIT NONE


    CONTAINS


    SUBROUTINE GenExcitations3(nI,iLut,tCountOnly,nSingles,nDoubles,Orbi,Orbj,Orba,Orbb,tAllExcitFound)
! This is the main routine called.  
! If tCountOnly is passed through as true, the single and double excitations are simply counted and passed back.
! If tCountOnly is false, the single and then double excitations in the form of Orbi,Orbj -> Orba,Orbb are
! found one by one.  
! The first single is found by passing in 0,0,0,0, and from then on, depending on the Orbi,Orbj,Orba and Orbb passed
! through, the next excitation in line will be found.
! When there are no more symmetry allowed excitations, tAllExcitFound becomes true.
        INTEGER :: nI(NEl),iLut(0:NIfD),Orbi,Orbj,Orba,Orbb,nSingles,nDoubles
        LOGICAL :: tCountOnly,tAllExcitFound
        
        IF(tCountOnly) THEN

            CALL CountExcitations(iLut,nI,nSingles,nDoubles)

        ELSEIF((Orbj.eq.0).and.(Orbb.eq.0)) THEN
! Generate singles, returning Orbi and Orba as non-zero, but keeping the others 0.        

            CALL GenSingleExcit(nI,iLut,Orbi,Orbj,Orba,Orbb)
            
            ! When the last single is input, generate the first double.

        ELSE
! This double is then passed in and the subsequent doubles found here.            

            CALL GenDoubleExcit(nI,iLut,Orbi,Orbj,Orba,Orbb,tAllExcitFound)

        ENDIF


    ENDSUBROUTINE GenExcitations3

    



    SUBROUTINE CountExcitations(iLut,nI,nSingleExcits,nDoubleExcits)
! This routine simply counts the excitations in terms of single and doubles from the nI determinant.    
        USE SymData, only: nSymLabels
        USE SystemData , only: ElecPairs
        USE GenRandSymExcitNUMod , only: PickElecPair,ConstructClassCounts 
        INTEGER :: nSingleExcits,nDoubleExcits,Symi,i,j,Spini,iLut(0:NIfD),nI(NEl)
        INTEGER :: iSpn,Elec1Ind,Elec2Ind,SymProduct
        INTEGER :: Syma,Symb,Spina,Spinb,StartSpin,EndSpin
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)


        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
! This sets up arrays containing the number of occupied and unoccupied in each symmetry.
! ClassCounts2(1,:)=No alpha occupied, ClassCounts2(2,:)=No Beta occupied.
! ClassCountsUnocc2(1,:)=No alpha unocc, ClassCounts2Unocc2(2,:)=No Beta unocc.
! The second index of these arrays referrs to the symmetry (0 -> 7).

! Only counting.  Run through each occupied orbital, and count the number of spin and symmetry allowed orbitals it 
! may be excited to.
        nSingleExcits=0
        nDoubleExcits=0
        do i=1,NEl
            Symi=INT(G1(nI(i))%Sym%S,4)
            IF((G1(nI(i))%Ms).eq.-1) Spini=2        ! G1(i)%Ms is -1 for beta, and 1 for alpha.
            IF((G1(nI(i))%Ms).eq.1) Spini=1         ! Translate this into 1 for alpha and 2 for beta.

! This electron in orbital of SymI and SpinI can only be excited to orbitals with the same spin and symmetry.                
! Then add in the number of unoccupied orbitals with the same spin and symmetry to which each electron may be excited.
            
            nSingleExcits=nSingleExcits+ClassCountUnocc2(Spini,Symi)

        enddo
! This is the end of the singles.        
        WRITE(6,*) 'Number of singles',nSingleExcits

! For the doubles, first pick an electron pair i,j.
! Based on these orbitals, run through each spin and each symmetry - take this to be orbital a.
! Multiply the number with these symmetries by the number of possible b orbitals which correspond.
! Do this for all a and then all i,j pairs.
        do i=1,ElecPairs
            CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,i)
! iSpn=2 for alpha beta pair, ispn=3 for alpha alpha pair and ispn=1 for beta beta pair.

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
                        ! There will exist a case where Orba = Orbb, want to remove this.
                        nDoubleExcits=nDoubleExcits+(ClassCountUnocc2(Spina,Syma)*(ClassCountUnocc2(Spinb,Symb)-1))
                    ELSE
                        nDoubleExcits=nDoubleExcits+(ClassCountUnocc2(Spina,Syma)*ClassCountUnocc2(Spinb,Symb))
                    ENDIF
                enddo

            enddo
        enddo

        WRITE(6,*) 'Number of doubles',nDoubleExcits

    ENDSUBROUTINE CountExcitations




    SUBROUTINE GenSingleExcit(nI,iLut,Orbi,Orbj,Orba,Orbb)
! Despite being fed four indices, this routine finds single excitations.  Orbi -> Orba. (Orbj and Orbb remain 0).
! Feeding in 0 indices indicates it is the first excitation that needs to be found.
! The single excitation goes from orbital i to a, from determinant nI to nJ.
! When the last single is found it then finds the first double excitation.
        USE SymData, only: nSymLabels
        INTEGER :: i,a,nI(NEl),Orbi,Orba,Symi,Finala,iLut(0:NIfD)
        INTEGER :: Orbj,Orbb,NoOcc
        LOGICAL :: tInitOrbsFound,tAllExcitFound
        INTEGER , SAVE :: OrbiIndex,OrbaIndex,Spini

        WRITE(6,*) 'Original Determinant',nI

        tInitOrbsFound=.false.

        IF((Orbi.eq.0).or.(Orba.eq.0)) THEN           ! Want to find the first excitation.

            OrbiIndex=1
            Orbi=nI(OrbiIndex)                              ! Take the first occupied orbital
            Symi=INT(G1(Orbi)%Sym%S,4)                      ! and find its spin and spat symmetries.
            IF((G1(Orbi)%Ms).eq.-1) Spini=2  
            IF((G1(Orbi)%Ms).eq.1) Spini=1  
            OrbaIndex=SymLabelCounts2(Spini,1,Symi+1)       ! Start considering a at the first allowed symmetry.

        ELSE
            Orbi=nI(OrbiIndex)                              ! Begin by using the same i as last time - check if there are any 
                                                            ! more possible excitations from this.

! At this stage, OrbaIndex is the a from the previous excitation.
            IF((OrbaIndex.eq.(nBasis/2)).or.(INT(G1(SymLabelList2(Spini,OrbaIndex+1))%Sym%S,4).ne.INT(G1(Orbi)%Sym%S,4))) THEN
! Either we're got to the final spin symmetry, or the next orbital after Orba does not have the same symmetry as Orbi.                
! Need to move onto the next i, and find a new a to match.
                OrbiIndex=OrbiIndex+1
                IF(OrbiIndex.le.NEl) THEN
                    Orbi=nI(OrbiIndex)
                    Symi=INT(G1(Orbi)%Sym%S,4)                  
                    IF((G1(Orbi)%Ms).eq.-1) Spini=2  
                    IF((G1(Orbi)%Ms).eq.1) Spini=1  
                    OrbaIndex=SymLabelCounts2(Spini,1,Symi+1)
                ENDIF

            ELSE
! There are more possible excitations from orbital a, simply check the next orbital after the current a.
                OrbaIndex=OrbaIndex+1
                Symi=INT(G1(Orbi)%Sym%S,4)           
            ENDIF
        ENDIF

        do while (.not.tInitOrbsFound)

            IF(OrbiIndex.gt.NEl) THEN
! If we've read in the last single, set orbi, orbj, orba, and orbb to 0 and call gendoubleexcit.        
                Orbi=0
                Orbj=0
                Orba=0
                Orbb=0
                CALL GenDoubleExcit(nI,iLut,Orbi,Orbj,Orba,Orbb,tAllExcitFound)
                EXIT
            ENDIF

! To find Orba, take the first in SymLabelList2 with the same symmetry and spin.                
! SymLabelCounts2(spin,1,symmetry) gives the index in SymLabelList2 where that spin and symmetry starts.                
            Orba=SymLabelList2(Spini,OrbaIndex)

! Need to also make sure orbital a is unoccupied, so make sure the orbital is not in nI.
            NoOcc=0
            do while (BTEST(iLut((Orba-1)/32),MOD((Orba-1),32))) 
! While this is true, Orba is occupied, so keep incrementing Orba until it is not.                    
                NoOcc=NoOcc+1
                Orba=SymLabelList2(Spini,OrbaIndex+NoOcc)
                IF(NoOcc.gt.SymLabelCounts2(Spini,2,Symi+1)) EXIT
            enddo

! Then check we have not overrun the symmetry block while skipping the occupied orbitals.                
            IF(INT(G1(Orba)%Sym%S,4).eq.Symi) THEN
! If not, then these are the new Orbi and Orba.                
                tInitOrbsFound=.true.
                OrbaIndex=OrbaIndex+NoOcc
            ELSE

! If we have, move onto the next occupied orbital i, no symmetry allowed single excitations exist from the first.                
                OrbiIndex=OrbiIndex+1
                IF(OrbiIndex.le.NEl) THEN
                    Orbi=nI(OrbiIndex)
                    Symi=INT(G1(Orbi)%Sym%S,4)                      ! and find its spin and spat symmetries.
                    IF((G1(Orbi)%Ms).eq.-1) Spini=2  
                    IF((G1(Orbi)%Ms).eq.1) Spini=1  
                    OrbaIndex=SymLabelCounts2(Spini,1,Symi+1)
                ENDIF
            ENDIF

        enddo
        WRITE(6,*) 'Excitation is from :',Orbi,' to ',Orba
        CALL FLUSH(6)


    ENDSUBROUTINE GenSingleExcit





    SUBROUTINE GenDoubleExcit(nI,iLut,Orbi,Orbj,Orba,Orbb,tAllExcitFound)
! This generates one by one, all possible double excitations.
! This involves a way of ordering the electron pairs i,j and a,b so that given an i,j and a,b we can find the next.
! The overall symmetry must also be maintained - i.e. if i and j are alpha and beta, a and b must be alpha and beta
! or vice versa.
        USE SystemData , only: ElecPairs
        USE GenRandSymExcitNUMod , only: PickElecPair 
        INTEGER :: nI(NEl),iLut(0:NIfD),Orbj,Orbi,Orba,Orbb,OrbbSpin,Syma,Symb
        INTEGER :: Elec1Ind,Elec2Ind,SymProduct,iSpn,Spinb
        INTEGER , SAVE :: ijInd,OrbaIndex,OrbbIndex,Spina
        LOGICAL :: tDoubleExcitFound,tFirsta,tFirstb,tNewij,tNewa,tAllExcitFound

        tDoubleExcitFound=.false.
        tFirsta=.false.
        tFirstb=.false.

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
            CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,ijInd)

            tNewij=.false.
! This becomes true when we can no longer find an allowed orbital a for this ij pair and we need to move onto the next.
            do while ((.not.tNewij).and.(.not.tDoubleExcitFound))                  ! This loop runs through the allowed a orbitals
                                                                                   ! until a double excitation is found.


                IF(tFirsta) THEN                    
! If this is the first double we are picking with this ij, we start with the alpha spin, unless i and j are both beta.
! There is no restriction on the symmetry for orbital a - although clearly the symmetry we pick determins b.
                    IF((iSpn.eq.2).or.(iSpn.eq.3)) Spina=1 
                    IF(iSpn.eq.1) Spina=2
                    OrbaIndex=1
                ENDIF

! If it is not the first, we have stored the previous spina and orba index - need to start with these and see                
! if any more double remain.
                Orba=SymLabelList2(Spina,OrbaIndex)

! The orbital chosen must be unoccupied.  This is just a test to make sure this is the case.
                do while (BTEST(iLut((Orba-1)/32),MOD((Orba-1),32))) 

! If not, we move onto the next orbital.                    
                    OrbaIndex=OrbaIndex+1

                    IF(OrbaIndex.gt.(nBasis/2)) THEN
! Have got to the end of that spin state, if ispn=2 the a orbital can be either alpha or beta, so we can just
! move onto the beta.
                        IF((iSpn.eq.2).and.(Spina.eq.1)) THEN
                            Spina=2
                            OrbaIndex=1
                            Orba=SymLabelList2(Spina,OrbaIndex)
                        ELSE
! We have got to the end of the possible a's - need to choose a new ij to start again.                            
                            tNewij=.true.
                            EXIT
                        ENDIF
                    ENDIF

! Otherwise the new orbital a is the first unoccupied orbital of allowed symmetry etc.                    
                    Orba=SymLabelList2(Spina,OrbaIndex)
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
                    Syma=INT(G1(Orba)%Sym%S,4)
                    Symb=IEOR(Syma,SymProduct)

! If this is the first time we've picked an orbital b for these i,j and a, begin at the start of the symmetry block.
! Otherwise pick up where we left off last time.
                    IF(tFirstb) THEN
                        OrbbIndex=SymLabelCounts2(Spinb,1,Symb+1)
                    ELSE
                        OrbbIndex=OrbbIndex+1
                    ENDIF

! If the new b orbital is still within the limits, check it is unoccupied and move onto the next orbital if it is.                    
                    IF((OrbbIndex.le.(nBasis/2)).and.(INT(G1(SymLabelList2(Spinb,OrbbIndex))%Sym%S,4).eq.Symb)) THEN
                        Orbb=SymLabelList2(Spinb,OrbbIndex)
                        do while ((BTEST(iLut((Orbb-1)/32),MOD((Orbb-1),32))).or.(Orbb.eq.Orba))
                            OrbbIndex=OrbbIndex+1
                            IF((OrbbIndex.gt.(nBasis/2)).or.(INT(G1(SymLabelList2(Spinb,OrbbIndex))%Sym%S,4).ne.Symb)) THEN
                                tNewa=.true.
                                EXIT
                            ENDIF
                            Orbb=SymLabelList2(Spinb,OrbbIndex)
                        enddo
                    ELSE
! If we have already gone beyond the symmetry limits by choosing the next b orbital, pick a new a orbital.                        
                        tNewa=.true.
                    ENDIF

! If we are moving onto the next a orbital, check we don't also need a new ij pair.                    
                IF(tNewa) THEN
                        OrbaIndex=OrbaIndex+1
                        tFirstb=.true.
                        IF(OrbaIndex.gt.(nBasis/2)) THEN
                            IF((iSpn.eq.2).and.(Spina.eq.1)) THEN
                                Spina=2
                                OrbaIndex=1
                                EXIT
                            ELSE
                                tNewij=.true.
                                ijInd=ijInd+1
                                IF(ijInd.gt.ElecPairs) THEN
                                    tAllExcitFound=.true.
                                ELSE
                                    EXIT
                                ENDIF
                            ENDIF
                        ELSE
                            EXIT
                        ENDIF
                    ENDIF

! If we get to here without exiting the loop, a set of four orbitals have been found. 
! Can exit the subroutine now.
                    tDoubleExcitFound=.true.

                enddo

            enddo

! This is the loop for new ij pairs - if we are choosing a new ij we are automatically choosing a new a and b also.            
            
            tFirsta=.true.
            tFirstb=.true.

        enddo

        Orbi=nI(Elec1Ind)
        Orbj=nI(Elec2Ind)

        WRITE(6,*) 'Excitation from : ',Orbi,Orbj,' to ',Orba,Orbb
        WRITE(6,*) 'These have symmetries : ',INT(G1(Orbi)%Sym%S,4),INT(G1(Orbj)%Sym%S,4),' to ',INT(G1(Orba)%Sym%S,4),INT(G1(Orbb)%Sym%S,4)
        CALL FLUSH(6)


    ENDSUBROUTINE GenDoubleExcit


END MODULE SymExcit3

