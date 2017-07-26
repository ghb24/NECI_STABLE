#include "macros.h"

module back_spawn_excit_gen

    use constants, only: dp, n_int, EPS, bits_n_int
    use SystemData, only: nel, G1, nbasis, tHPHF, NMAXX, NMAXY, NMAXZ, &
                          tOrbECutoff, OrbECutoff, nOccBeta, nOccAlpha
    use bit_rep_data, only: niftot
    use SymExcitDataMod, only: excit_gen_store_type, SpinOrbSymLabel, kPointToBasisFn
    use bit_reps, only: test_flag, get_initiator_flag
    use FciMCData, only: pSingles, projedet, pDoubles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use excit_gens_int_weighted, only: gen_single_4ind_ex, select_orb_sing, &
                                pick_weighted_elecs, get_paired_cc_ind, select_orb, &
                                pgen_select_orb, pgen_weighted_elecs, pgen_single_4ind
    use excit_gen_5, only: gen_double_4ind_ex2, pick_a_orb, pgen_select_a_orb, &
                        calc_pgen_4ind_weighted2
    use CalcData, only: t_back_spawn_flex, t_back_spawn_occ_virt, t_back_spawn, &
                        occ_virt_level, t_back_spawn_flex_option, t_back_spawn_option
    use GenRandSymExcitNUMod, only: ClassCountInd, RandExcitSymLabelProd, &
                        CreateExcitLattice
    use back_spawn, only: check_electron_location, pick_virtual_electrons_double, & 
                          pick_occupied_orbital_single, pick_virtual_electron_single, &
                          pick_occupied_orbital, pick_second_occupied_orbital
    use get_excit, only: make_single, make_double
    use Determinants, only: write_det, get_helement

!     use hphf_integrals, only: hphf_off_diag_helement

#ifdef __DEBUG 
    use SystemData, only: tNoFailAb
#endif

    implicit none

contains

    ! to disentangle the mess of excitation generators also implement an new
    ! one for the UEG and the hubbard model seperately 
    subroutine gen_excit_back_spawn_ueg(nI, ilutI, nJ, ilutJ, exFlag, ic, &
            ExcitMat, tParity, pgen, HelGen, store, run) 
        ! if back-spawn and ueg is turned on point to this excitation 
        ! generator! check if we hit all the relevant parts in the code though
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store 
        integer, intent(in), optional :: run
        character(*), parameter :: this_routine = "gen_excit_back_spawn_ueg"
 
        logical :: temp_back_spawn
        ! so now pack everything necessary for a ueg excitation generator. 

        ilutJ(0) = -1 
        HElGen = 0.0_dp 
        ic = 2

        ! this function gets pointed to if tUEG, t_back_spawn_flex and tLatticeGens
        ! is set 
        ! BUT: again: we also need to take care of the HPHF keyword.. 
        ! which is a mess 
        temp_back_spawn = (t_back_spawn_flex .and. .not. &
            test_flag(ilutI, get_initiator_flag(run))) 

        ! implement it for now without this tNoFailAb flag
        ASSERT(.not. tNoFailAb)

        if (temp_back_spawn) then 

            call gen_double_back_spawn_ueg(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, & 
                pgen) 

        else 
            ! do i want to rewrite the old on or just reuse? for now reuse: 
            call CreateExcitLattice(nI, ilutI, nJ, tParity, ExcitMat, pgen, run) 

        end if



    end subroutine gen_excit_back_spawn_ueg

    subroutine gen_double_back_spawn_ueg(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, & 
                pgen, run) 
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
        integer, optional :: run
        character(*), parameter :: this_routine = "gen_double_back_spawn_ueg"

        integer :: kaxrange, kayrange, kazrange, ielecinexcitrange, Elec1, & 
                   Elec2, Elec1Ind, Elec2Ind, iSpn, hole2basisnum, elecs(2), &
                   loc, temp_run, ki(3), kj(3), ka(3), kb(3), kb_ms, TestEnergyB, &
                   iSpinIndex, ChosenUnocc, hole1basisnum
        real(dp) :: r(2), x, pAIJ, dummy
        logical :: tAllowedExcit

        if (present(run)) then 
            temp_run = run 
        else 
            temp_run = 1
        end if

        kaxrange = 0
        kayrange = 0
        kazrange = 0
        ielecinexcitrange = 0

        ! i know here that back-spawn is on and this is a non-inititator: 
        ! so for the beginning implementation we pick the two electrons 
        ! randomly 
        ! i should make this to a function below: maybe with the hubbard 
        ! flag as additional input to provide us with two independent 
        ! electrons in any order possible!
        DO
            r(1) = genrand_real2_dSFMT()
            Elec1=INT(r(1)*NEl+1)
            DO
                r(2) = genrand_real2_dSFMT()
                Elec2=INT(r(2)*NEl+1)
                IF(Elec2.ne.Elec1) EXIT
            ENDDO
            Elec1Ind=Elec1
            Elec2Ind=Elec2
            IF((G1(nI(Elec1Ind))%Ms.eq.-1).and.(G1(nI(Elec2Ind))%Ms.eq.-1)) THEN
                iSpn=1
            ELSE
                IF((G1(nI(Elec1Ind))%Ms.eq.1).and.(G1(nI(Elec2Ind))%Ms.eq.1)) THEN 
                    iSpn=3
                ELSE
                    iSpn=2
                ENDIF
            ENDIF
        ENDDO

        ! i should make it more efficient to pick the hole.. 
        hole2basisnum = 0
        elecs = [Elec1Ind, Elec2Ind]

        loc = check_electron_location(elecs, 2, temp_run)
        
        if (loc == 2) then 
            do 
                x = genrand_real2_dSFMT() 
                if (iSpn == 2) then 
                    ChosenUnocc = int(nBasis * x) + 1

                else 
                    ChosenUnocc=2*(INT(nBasis/2*x)+1) & ! 2*(a number between 1 and nbasis/2) gives the alpha spin
                    & -(1-(iSpn/3)) ! alpha numbered even, iSpn/3 returns 1 for alpha/alpha, 0 for beta/beta
                end if
                
                if (IsNotOcc(ilutI, ChosenUnocc)) exit 
            end do
            pAIJ = 1.0_dp / real(nBasis - nel, dp) 

        else 
            ! i think i need a new one for the ueg also.. since there is not 
            ! such a spin restriction as in the hubbard model
            ! i think just using the "normal" occupied picker should be fine
            call pick_occupied_orbital_ueg(nI, elecs, ispn, temp_run, pAIJ, &
                     dummy, ChosenUnocc)

        end if

        ! i guess this is not necessary, but anyway..
        Hole1BasisNum=ChosenUnocc
        IF((Hole1BasisNum.le.0).or.(Hole1BasisNum.gt.nBasis)) THEN
            CALL Stop_All("CreateDoubExcitLattice","Incorrect basis function generated") 
        ENDIF

        ! kb is now uniquely defined
        ki=G1(nI(Elec1Ind))%k
        kj=G1(nI(Elec2Ind))%k
        ka=G1(Hole1BasisNum)%k
        kb=ki+kj-ka

        ! Find the spin of b
        IF(iSpn.eq.2)THEN ! alpha/beta required, therefore b has to be opposite spin to a
            kb_ms=1-((G1(Hole1BasisNum)%Ms+1)/2)*2
        ELSE ! b is the same spin as a
            kb_ms=G1(Hole1BasisNum)%Ms
        ENDIF

        ! Is kb allowed by the size of the space?
        ! Currently only applies when NMAXX etc. are set by the CELL keyword
        ! Not sure what happens when an energy cutoff is set
        tAllowedExcit=.true.
        IF(ABS(kb(1)).gt.NMAXX) tAllowedExcit=.false.
        IF(ABS(kb(2)).gt.NMAXY) tAllowedExcit=.false.
        IF(ABS(kb(3)).gt.NMAXZ) tAllowedExcit=.false.
        TestEnergyB=kb(1)**2+kb(2)**2+kb(3)**2
        IF(tOrbECutoff.and.(TestEnergyB.gt.OrbECutoff)) tAllowedExcit=.false.
        IF(.not.tAllowedExcit) THEN
            nJ(1)=0
            RETURN
        ENDIF
               
        iSpinIndex=(kb_ms+1)/2+1
        Hole2BasisNum=kPointToBasisFn(kb(1),kb(2),kb(3),iSpinIndex)
       
        IF(Hole2BasisNum==-1.or.Hole1BasisNum.eq.Hole2BasisNum) THEN
            nJ(1)=0 
            RETURN
        ENDIF

        ! Is b occupied?
        IF(BTEST(iLutI((Hole2BasisNum-1)/bits_n_int),MOD(Hole2BasisNum-1,bits_n_int))) THEN
        !Orbital is in nI. Reject.
            tAllowedExcit=.false.
        ENDIF
        
        IF(.not.tAllowedExcit) THEN
            nJ(1)=0
            RETURN
        ENDIF

        ! Find the new determinant
        call make_double (nI, nJ, elec1ind, elec2ind, Hole1BasisNum, &
                          Hole2BasisNum, ExcitMat, tParity)
                
        ! we knoe it is back-spawn + UEG remember! so pAIJ is already above
        !Calculate generation probabilities
        if (loc /= 2) then
            IF (iSpn.eq.2) THEN
                ! otherwise take the pAIJ set above
                pAIJ=1.0_dp/(nBasis-Nel)
            ELSEIF (iSpn.eq.1) THEN
                pAIJ=1.0_dp/(nBasis/2-nOccBeta)
            ELSE
                !iSpn = 3
                pAIJ=1.0_dp/(nBasis/2-nOccAlpha)
            ENDIF
        end if
        ! Note, p(b|ij)=p(a|ij) for this system
        pGen=2.0_dp/(NEl*(NEl-1))*2.0_dp*pAIJ

        ! pgens at the end.. 

    end subroutine gen_double_back_spawn_ueg

    
    ! also write a wrapper-like routine for an excitation generator if 
    ! back-spawn is activated.. to not mess up all the old functions too much. 
    subroutine gen_excit_back_spawn(nI, ilutI, nJ, ilutJ, exFlag, ic, &
            ExcitMat, tParity, pgen, HelGen, store, run) 
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store 
        integer, intent(in), optional :: run
        character(*), parameter :: this_routine = "gen_excit_back_spawn"

        logical :: temp_back_spawn 
#ifdef __DEBUG
        HElement_t(dp) :: temp_hel
        real(dp) :: pgen2
#endif
        ! check the non-initiator criteria beforehand 
        ! i also have to consider that back-spawn gets turned on later on 
        ! so i have to check if back-spawn is active already or not..
        temp_back_spawn = ((t_back_spawn_flex .or. t_back_spawn) .and. .not. &
            test_flag(ilutI, get_initiator_flag(run)))

        if (temp_back_spawn) then 

            ! otherwise use custom made ones
            if (genrand_real2_dSFMT() < pSingles) then 
                
                ic = 1
                call gen_single_back_spawn(nI, ilutI, run, nJ, ilutJ, ExcitMat, &
                    tParity, pgen)
                pgen = pgen * pSingles 

            else 
                
                ic = 2 
                call gen_double_back_spawn(nI, ilutI, run, nJ, ilutJ, ExcitMat, &
                    tParity, pgen) 
                pgen = pgen * pDoubles 

            end if

! #ifdef __DEBUG
!             if (.not. IsNullDet(nJ)) then
!                 pgen2 = calc_pgen_back_spawn(nI, ilutI, ExcitMat, ic, run) 
!                 if (abs(pgen - pgen2) > 1.0e-6_dp) then 
!                     if (tHPHF) then 
!                         print *, "due to circular dependence, no matrix element calc possible!"
! !                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
!                         temp_hel = 0.0_dp
!                     else
!                         temp_hel = get_helement(nI,nJ,ilutI,ilutJ)
!                     endif
! 
!                     write(6,*) 'Calculated and actual pgens differ. for non-initiator'
!                     write(6,*) 'This will break HPHF calculations'
!                     call write_det(6, nI, .false.)
!                     write(6, '(" --> ")', advance='no')
!                     call write_det(6, nJ, .true.)
!                     write(6,*) 'Excitation matrix: ', ExcitMat(1,1:ic), '-->', &
!                                ExcitMat(2,1:ic)
!                     write(6,*) 'Generated pGen:  ', pgen
!                     write(6,*) 'Calculated pGen: ', pgen2
!                     write(6,*) 'matrix element: ', temp_hel
!                     call stop_all(this_routine, "Invalid pGen")
!                 end if
!             end if
! #endif
        else 
 
            ! do the "normal" excitation type if it is an initiator
            if (genrand_real2_dSFMT() < pSingles) then 

                ic = 1 
                call gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, &
                                     tParity, pGen)
                pgen = pgen * pSingles
                
            else 

                ic = 2
                call gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ExcitMat, &
                    tParity, pGen)
                pgen = pgen * pDoubles

            end if
! #ifdef __DEBUG
!             if (.not. IsNullDet(nJ)) then
!                  pgen2 = calc_pgen_4ind_weighted2(nI, ilutI, ExcitMat, ic)
!                 if (abs(pgen - pgen2) > 1.0e-6_dp) then
!                     if (tHPHF) then
!                         print *, "due to circular dependence, no matrix element calc possible!"
! !                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
!                         temp_hel = 0.0_dp
!                     else
!                         temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
!                     end if
! 
!                     write(6,*) 'Calculated and actual pgens differ. initiator'
!                     write(6,*) 'This will break HPHF calculations'
!                     call write_det(6, nI, .false.)
!                     write(6, '(" --> ")', advance='no')
!                     call write_det(6, nJ, .true.)
!                     write(6,*) 'Excitation matrix: ', ExcitMat(1,1:ic), '-->', &
!                                ExcitMat(2,1:ic)
!                     write(6,*) 'Generated pGen:  ', pgen
!                     write(6,*) 'Calculated pGen: ', pgen2
!                     write(6,*) 'matrix element: ', temp_hel
!                     call stop_all(this_routine, "Invalid pGen")
!                 end if
!             end if
! #endif
       end if

    end subroutine gen_excit_back_spawn

    subroutine gen_single_back_spawn(nI, ilutI, run, nJ, ilutJ, ex, tPar, pgen) 
        ! specialised single excitation routine for the back-spawn method
        ! for the moment i still have to decide, which back-spawn method is 
        ! in use.. 
        integer, intent(in) :: nI(nel) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: run
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_single_back_spawn"

        integer :: elec, src, cc_index, loc, tgt
        real(dp) :: pgen_elec

        ! depending on the method we pick electrons accordingly
        if (t_back_spawn_flex) then 
            elec = 1 + floor(genrand_real2_dSFMT() * nel)

            loc = check_electron_location([nI(elec),0], 1, run)

            pgen_elec = 1.0_dp/real(nel, dp)
        else
            call pick_virtual_electron_single(nI, run, elec, pgen_elec)
        end if

        src = nI(elec)

        cc_index = ClassCountInd (get_spin(src), SpinOrbSymLabel(src), &
                                  G1(src)%Ml)


        ! i have to make the logic easier here at some point.. 
        ! we just need to do some extensive testing and decide on one 
        ! back-spawn method and abandon the rest.. 
        ! i hope this logic is correct: otherwise check on the bottom
        if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. (&
            (loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2))) then 

            call pick_occupied_orbital_single(nI, ilutI, src, cc_index, run, pgen, tgt)

        else 

            tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 

        end if
! 
!         if (t_back_spawn_occ_virt) then 
!             call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)
! 
!         else if (t_back_spawn_flex) then 
! 
!             if (loc == 2 .and. occ_virt_level /= -1) then 
!                 call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 
! 
!             else 
!                 if (occ_virt_level == 2) then 
!                     call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 
!                 else
!                     tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
!                 end if 
!             end if
!         else
!             tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
!         end if

        if (tgt == 0) then 
            nJ(1) = 0
            pgen = 0.0_dp 
            return
        end if

        call make_single(nI, nJ, elec, tgt, ex, tPar)

        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)

        pgen = pgen * pgen_elec

    end subroutine gen_single_back_spawn

    subroutine gen_double_back_spawn(nI, ilutI, run, nJ, ilutJ, ex, tPar, pgen) 
        ! the double excitation routine for the back-spawn method 
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: run
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_double_back_spawn"

        integer :: elecs(2), sym_product, ispn, sum_ml, src(2), loc, cc_a, cc_b,&
                   orbs(2)
        real(dp) :: int_cpt(2), cum_sum(2), sum_pair(2), cpt_pair(2), &
                    cum_arr(nbasis)

        if (t_back_spawn_flex) then 
            ! Pick the electrons in a weighted fashion
            call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, &
                                 pgen)

            loc = check_electron_location(src, 2, run)

        else
            call pick_virtual_electrons_double(nI, run, elecs, src, ispn,&
                                                sum_ml, pgen)

            ! here it could be that there are no valid electrons.. 
            if (elecs(1) == 0) then 
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            ! TODO!! thats the problem here! there are circular dependencies 
            ! when I call that from here.. so i guess i need to setup the 
            ! back_spawn excitation routines in a seperate file 
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                             SpinOrbSymLabel(src(2)))

            ! indicate no flex option is used
            loc = -1

        end if

        ! at some point i have to make this whole logic easier here.. 
        ! in the end i think i have to remove some possibilities and stick to 
        ! one back-spawn method. 
        if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. (&
            (loc == 1 .and. occ_virt_level /= -1) .or. (loc == 2) .or. &
            (loc == 0 .and. occ_virt_level >= 1)))) then 

            call pick_occupied_orbital(nI, src, ispn, run, int_cpt(1), cum_sum(1), &
                                        orbs(1))

        else 

            orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
            
        end if

! 
!         if (t_back_spawn_temp ) then
!  
!             ! if we have one of the electrons in the occupied manifold 
!             ! pick atleast one hole also from this manifold to not increase 
!             ! the excitation level!
!            
!             call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
! 
!         else if (t_back_spawn_flex ) then
!             ! now we have to decide on the flex-spawn + occ-virt implo:
!             if (loc == 1) then 
!                 ! the new option to excite on level
!                 if (occ_virt_level == -1) then
!                     orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!                 else
!                     call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!                 end if
! 
! 
!             else if (loc == 2) then 
!                 ! then we always pick an occupied orbital
!                 call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!             else 
!                 ! depending on the occ_virt_level
!                 if (occ_virt_level < 1) then 
!                     ! just pick normal:
!                     orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!                 else 
!                     ! otherwise if it is 1 or 2, we want to pick the (a) also
!                     ! from the occupied manifold
! 
!                     call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!                 end if
!             end if
!         else 
!             orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!         end if


        if (orbs(1) /= 0) then
            cc_a = ClasSCountInd(orbs(1))
            cc_b = get_paired_cc_ind(cc_a, sym_product, sum_ml, iSpn)

            if (t_back_spawn_flex .and.((loc == 2 .and. occ_virt_level /= -1) .or. &
                (occ_virt_level == 2))) then 

                call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
                    run, int_cpt(2), cum_sum(2), orbs(2))

            else 

                orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
                                  cum_sum(2))
            end if

            ASSERT((.not. (is_beta(orbs(2)) .and. .not. is_beta(orbs(1)))))

!             if (t_back_spawn_flex  .and. .not. t_temp_init) then 
!                 ! in this case i have to pick the second orbital also from the 
!                 ! occupied list, but now also considering symmetries
!                 if (loc == 2) then 
!                     ! now also mix the occ-virt with this back-spawning
!                     ! for loc 2 always do it
!                     ! except specified by occ_virt_level= -1
!                     if (occ_virt_level == -1) then
!                         orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!                     else
!                         call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                             int_cpt(2), cum_sum(2), orbs(2))
!                     end if
! 
!                 else if (loc == 1 .and. occ_virt_level == 2) then 
!                     call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                         int_cpt(2), cum_sum(2), orbs(2))
! 
!                 else if (loc == 0 .and. occ_virt_level == 2) then 
!                     call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                         int_cpt(2), cum_sum(2), orbs(2))
!                 else
!                     orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!                 end if
! 
!             else
! 
!                 orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!             end if

        end if

        if (any(orbs == 0)) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! can i exit right away if this happens??
        ! i am pretty sure this means 
        if (any(cum_sum < EPS)) then 
           cum_sum = 1.0_dp
           int_cpt = 0.0_dp
        end if
        
        ! only on parallel excitations.. and symmetric exciation generator is 
        ! turned off for now in the back-spawning
        if (is_beta(orbs(1)) .eqv. is_beta(orbs(2)))  then
            if (t_back_spawn_occ_virt .or. (t_back_spawn_flex .and. (& 
                (loc == 1 .and. (occ_virt_level == 0 .or. occ_virt_level == 1)) &
                .or. (loc == 2 .and. occ_virt_level == -1) .or. &
                (loc == 0 .and. occ_virt_level == 1)))) then 

                if (any(orbs(2) == projedet(:,part_type_to_run(run)))) then
                   ! if (b) is also in the occupied manifold i could have 
                    ! picked the other way around.. 
                    ! with the same uniform probability: 
                    cpt_pair(1) = int_cpt(1)
                    sum_pair(1) = cum_sum(1) 
                    ! and then (a) would have been picked according to the 
                    ! "normal" procedure
                    call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                 cpt_pair(2), sum_pair(2))
                else
                    ! if (b) is not in the occupied this does not work or? 
                    ! since i am forcing (a) to be in the occupied.. 
                    ! so remove this pgen:
                    cpt_pair = 0.0_dp
                    sum_pair = 1.0_dp
                end if

            else if (t_back_spawn_flex .and. ((loc == 0 .and. occ_virt_level == 2) & 
                .or. (loc == 1 .and. occ_virt_level == 2) .or. & 
                (loc == 2 .and. occ_virt_level /= -1))) then 

                cpt_pair = int_cpt 
                sum_pair = cum_sum 
                

            else

                ! otherwise "normal"
                call pgen_select_a_orb(ilutI, src, orbs(2), iSpn, cpt_pair(1), &
                                       sum_pair(1), cum_arr, .false.)
                call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                     cpt_pair(2), sum_pair(2))

            end if 

        else 

            cpt_pair = 0.0_dp
            sum_Pair = 1.0_dp
                
        end if

        if (any(sum_pair < EPS)) then 
            cpt_pair = 0.0_dp
            sum_pair = 1.0_dp
        end if

        pgen = pgen * (product(int_cpt) / product(cum_sum) + &
                       product(cpt_pair) / product(sum_pair))

        ! And generate the actual excitation
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, tpar)
        ilutJ = ilutI
        clr_orb (ilutJ, src(1))
        clr_orb (ilutJ, src(2))
        set_orb (ilutJ, orbs(1))
        set_orb (ilutJ, orbs(2))


    end subroutine gen_double_back_spawn

    function calc_pgen_back_spawn_hubbard(nI, ilutI, ex, ic, run) result(pgen)
        integer, intent(in) :: nI(nel), ex(2,2), ic, run
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn"

        ! the hubbard pgen recalculation is pretty similar to the ueg one
        ! below..

    end function calc_pgen_back_spawn_hubbard

    function calc_pgen_back_spawn_ueg(nI, ilutI, ex, ic, run) result(pgen)
        integer, intent(in) :: nI(nel), ex(2,2), ic, run
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn"

        real(dp) :: cum_sum 
        integer :: dummy_orb, ispn, loc, src(2)

        ! i have to write a generation probability calculator for the hphf 
        ! implementation this should only be called if it is a non-init 
        ! and back-spawn is on.. otherwise just use the already provided 
        ! ones
        if (test_flag(ilutI, get_initiator_flag(run))) then 
            
            ! i have to do some symmetry setup beforehand.. 
            ! or i do it by hand to avoid the unnecessary overhead.. 
!             call calc_pgen_symrandexcit2(nI, ex, 2, 
        else 
            ! do the back-spawn pgen.. 
            ! elec were picked randomly(but in both orders!) 
            ! and then we have to check the electron location, to know if 
            ! there was a restriction on the first orbitals 
            loc = check_electron_location(src, 2, run) 

            ! and with the implementation right now it is only restricted 
            ! if both were in the occup
            if (loc == 2) then 
            
                ! i could just call  the orbital picking routine.. 
                ! although this is inefficient.. 
                ! and i have to see if it would have been possible the 
                ! other way around too.. 
                call pick_occupied_orbital(nI, src, ispn, run, pgen, cum_sum, &
                    dummy_orb)

            else 
                ! otherwise the orbital (a) is free 
                ! so it is the number of available orbitals, depending on 
                ! the spin 
                
            end if
        end if

    end function calc_pgen_back_spawn_ueg

    function calc_pgen_back_spawn(nI, ilutI, ex, ic, run) result(pgen)
        ! to use HPHF keyword and also the test if the pgens are correct 
        ! or just to be able and to be sure and more save i need a way to 
        ! recalculate the generation probability also for a back-spawn 
        ! excitation from a non-initiator determinant
        ! so although thats a hassle implement it, otherwise i cannot be 
        ! quite sure about the method
        integer, intent(in) :: nI(nel), ex(2,2), ic, run
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn"

        integer :: dummy, ssrc, stgt, cc_index, src(2), tgt(2), dummy_elecs(2), &
                   dummy_orbs(2), ispn, loc, sum_ml, sym_prod, cc_a, cc_b, &
                   dummy_ispn, dummy_sum_ml
        real(dp) :: elec_pgen, orb_pgen, int_cpt(2), cum_sum(2), cpt_pair(2), &
                    sum_pair(2), cum_arr(nbasis)
        logical :: t_gen_list, t_in_ref, t_par
        ! i should only call this function when i am on a non-inititator or?
        ! in the HPHF framework i could end up here also on a non-initiator
        ! so here i should then go to 4ind-2 pgen calculator if it is a 
        ! inititator 

        if (test_flag(ilutI, get_initiator_flag(run))) then 

            pgen = calc_pgen_4ind_weighted2(nI, ilutI, ex, ic)

        else 
            ! in the hphf framework we definetly only end up here if the 
            ! back-spawn method is already turned on.. but is this 
            ! ensured in the other places the function might get called? 
            ! decision: if we call this function we want to get the pgen in 
            ! the back-spawn method, so ensure outside if it is turned on 
            ! or not. 

            if (ic == 1) then 

                ! depending on the implementation
                ! here since i want to calculate the real pgen if back_spawn 
                ! is or would be on, i should work with the 
                ! _option keywords..
                ssrc = ex(1,1)
                stgt = ex(2,1)

                if (t_back_spawn_flex_option) then 
                    elec_pgen = 1.0_dp / real(nel, dp) 

                    loc = check_electron_location([ssrc,0], 1, run)


                else 
                    ! it is slow anyway.. so just do the dummy implementation
                    ! as Ali mentioned in the HPHF implementation it is not 
                    ! that common to have a doubly connected HPHF det

                    call pick_virtual_electron_single(nI, run, dummy, elec_pgen) 


                    loc = -1

                end if


                if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex_option .and.&
                    ((loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2))) then 


                    ! the symmetry of the orbital is known after all
                    cc_index = ClassCountInd(get_spin(stgt), SpinOrbSymLabel(stgt), &
                        G1(stgt)%Ml)

                    ! reuse the routine.. 
                    call pick_occupied_orbital_single(nI, ilutI, ssrc, cc_index, &
                        run, orb_pgen, dummy)

                else 


                    ! reuse the other pgen single function 
                    ! but i have to multiply out the electron pgen
                    orb_pgen = pgen_single_4ind(nI, ilutI, ssrc, stgt) * &
                        real(nel, dp)


                end if

                pgen = pSingles * elec_pgen * orb_pgen


            else if (ic == 2) then 


                src = ex(1,:)
                tgt = ex(2,:)

                ! get all the necessary symmetry info at the beginning
                if (is_beta(src(1)) .eqv. is_beta(src(2))) then
                    if (is_beta(src(1))) then
                        iSpn = 1
                    else
                        iSpn = 3
                    end if
                else
                    iSpn = 2
                end if
                sum_ml = sum(G1(src)%ml)

                sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                                 SpinOrbSymLabel(src(2)))

                if (t_back_spawn_flex_option) then 

                    elec_pgen = pgen_weighted_elecs(nI, src)

                    loc = check_electron_location(src, 2, run)


                else 

                    call pick_virtual_electrons_double(nI, run, dummy_elecs, &
                        dummy_orbs, dummy_ispn, dummy_sum_ml, elec_pgen) 


                    loc = -1

                end if

                ! i have to be careful with the orbitals.. 
                ! because i guess those get ordered.. and i have to check if 
                ! it is possible to have picked the orbitals in either 
                ! order.. 
                ! ok.. and there is the restriction to pick a beta orbital 
                ! in src(1) first if it is a anti-parallel excitation
                ! ok this means atleast src(1) is definetly the first picked 
                ! orbital.. is this also the case in HPHF?? argh i hate that 
                ! stuff 
                ! do this testing here once
                t_in_ref = (any(tgt(2) == projedet(:,part_type_to_run(run))))
                t_par = (is_beta(tgt(1)) .eqv. is_beta(tgt(2)))

                ! for some back_spawn_flex it can happen that we have no 
                ! restrictions on the orbitals.. but thats rare i guess..
                if (t_back_spawn_occ_virt .or. (t_back_spawn_flex_option .and. &
                    ((loc == 1 .and. occ_virt_level /= -1) .or. loc == 2 .or. &
                    (loc == 0.and. occ_virt_level >= 1)))) then 

                    call pick_occupied_orbital(nI, src, ispn, run, int_cpt(1),&
                        cum_sum(1), dummy_orbs(1))

                    ! i can atleast do some stuff for picking it the other
                    ! way or?
                    ! because if it is a parallel spin-excitations and orbital
                    ! (b) is in the reference the probs p(b|ij) are the same 
                    if (t_par .and. t_in_ref) then
                        cpt_pair(1) = int_cpt(1) 
                        sum_pair(1) = cum_sum(1) 

                    else 
                        cpt_pair = 0.0_dp
                        sum_pair = 1.0_dp
                        ! maybe i could set a flag that it does not have to 
                        ! be recalced otherwise.. 
                    end if

                    ! here i should go on to calculate p(b|aij) since in the 
                    ! other case we have everything.. 

                    if (t_back_spawn_flex_option .and. (&
                        (loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2)) then

                        ! this is the case where orbital (b) is also restricted

                        cc_a = ClassCountInd(tgt(1)) 
                        cc_b = get_paired_cc_ind(cc_a, sym_prod, sum_ml, ispn)

                        call pick_second_occupied_orbital(nI, src, cc_b, tgt(1), &
                            ispn, run, int_cpt(2), cum_sum(2), dummy_orbs(2))

                        ! and ofc both probs are the same if the spin-fits
                        if (t_par) then 
                            cpt_pair = int_cpt 
                            sum_pair = cum_sum 
                        else
                            cpt_pair = 0.0_dp
                            sum_pair = 1.0_dp
                        end if

                    else

                        call pgen_select_orb(ilutI, src, tgt(1), tgt(2), int_cpt(2), &
                            cum_sum(2))

                        ! in this case (a) was restricted but (b) was not 
                        ! so check if the other way would have been 
                        ! possible 
                        if (t_par .and. t_in_ref) then 
                            ! p(b|ij) already set above 
                            call pgen_select_orb(ilutI, src, tgt(2), tgt(1), &
                                cpt_pair(2), sum_pair(2)) 

                        else
                            cpt_pair = 0.0_dp
                            sum_pair = 1.0_dp
                        end if

                    end if

                else 

                    call pgen_select_a_orb(ilutI, src, tgt(1), ispn, int_cpt(1),&
                        cum_sum(1), cum_arr, .true.) 

                    ! but i guess in this case i can already cover all the 
                    ! pgen.. 
                    if (int_cpt(1) > EPS) then 
                        call pgen_select_orb(ilutI, src, tgt(1), tgt(2), &
                            int_cpt(2), cum_sum(2))

                        t_gen_list = .false.
                    else
                        t_gen_list = .false. 
                        int_cpt = 0.0_dp
                        cum_sum = 1.0_dp
                    end if
                    ! otherwise i can pick orbital (a) freely.. which also 
                    ! means that i definetly can pick (b) freely and do it 
                    ! the other way around.. 

                    call pgen_select_a_orb(ilutI, src, tgt(2), ispn, cpt_pair(1),&
                        sum_pair(1), cum_arr, t_gen_list)

                    if (cpt_pair(1) > EPS) then 
                        call pgen_select_orb(ilutI, src, tgt(2), tgt(1), &
                            cpt_pair(2), sum_pair(2))
                    else 
                        cpt_pair = 0.0_dp
                        sum_pair = 1.0_dp
                    end if
                end if 

                ! how do we handle this now?? 
                ! i only want to handle these cases, which i have not 
                ! recalculated above already. especially for the p(b|ij) 
                ! probability.. 

                ! now i have to figure p(ab), p(ba) stuff.. and implement this 
                ! above.. annoying..



            else 

                pgen = 0.0_dp

            end if
        end if
        

    end function calc_pgen_back_spawn

end module back_spawn_excit_gen
