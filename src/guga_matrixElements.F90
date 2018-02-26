#include "macros.h"
! GUGA module containg as much matrix element calculation functionality as 
! possible. 
#ifndef __CMPLX
module guga_matrixElements
    ! used modules: 
    use SystemData, only: nEl, nBasis, ECore
    use constants, only: dp, n_int, hel_zero
    use bit_reps, only: niftot, decode_bit_det, nifd
    use OneEInts, only: GetTMatEl
    use Integrals_neci, only: get_umat_el
    use guga_bitRepOps, only: isDouble, calcB_vector_nI, isProperCSF_nI, &
                            extract_matrix_element
    use util_mod, only: binary_search
    use guga_data, only: projE_replica
    use bit_rep_data, only: nifdbo
    use ParallelHelper, only: iprocindex

    ! variable declarations:
    implicit none
    
    ! interfaces
    !interface calcDiagMatEleGUGA
    !    ! general diagonal matrix element calculations:
    !    ! chosen depending on input
    !    module procedure calcDiagMatEleGUGA_nI
    !    module procedure calcDiagMatEleGuga_ilut
    !end interface calcDiagMatEleGUGA

    !interface calcDiagExchangeGUGA
    !    ! general diagonal exchange contribution calculation:
    !    module procedure calcDiagExchangeGUGA_nI
    !    module procedure calcDiagExchangeGUGA_ilut
    !end interface calcDiagExchangeGUGA

!     interface calc_off_diag_guga
!         module procedure calc_off_diag_guga_ref
!         module procedure calc_off_diag_guga_gen
!     end interface calc_off_diag_guga

contains

    function calc_off_diag_guga_ref_list(ilut, run, exlevel) result(hel)
        ! calculated the off-diagonal element connected to the reference
        ! determinant only. 
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in), optional :: run
        integer, intent(out), optional :: exlevel
        HElement_t(dp) :: hel
        character(*), parameter :: this_routine = "calc_off_diag_guga_ref"

        integer :: pos, ind, nExcit
        ! have the list of conncected dets to ilutRef stored persistently 
        ! so only need to search if ilut is in this list and return 
        ! the corresponing matrix element 
        ! for now change that to the first replica reference...
        ! since this is done also in the fcimc_helper in the sumEcontrib...
        ! im not quite sure why this is done in such a way in the main routine.
        if (present(run)) then
            ind = run
        else
            ind = 1
        end if

        ASSERT(allocated(projE_replica(ind)%projE_ilut_list))

        nExcit = projE_replica(ind)%num_entries

        pos = binary_search(projE_replica(ind)%projE_ilut_list(0:nifd,1:nExcit), &
            ilut(0:nifd))

        if (pos > 0) then
            ! if found output the matrix element 
            hel = projE_replica(ind)%projE_hel_list(pos)
            if (present(exlevel)) then
                exlevel = projE_replica(ind)%exlevel(pos)
            end if
        else 
            ! otherwise its zero
            hel = hel_zero
            if (present(exlevel)) then
                ! which value should i give exlevel in this case? 3,-1 ..
                ! have to deal with it outside..
                exlevel = -1
            end if
        end if

    end function calc_off_diag_guga_ref_list


    function calcDiagMatEleGuga_nI(nI) result(hel_ret)
        ! calculates the diagonal Hamiltonian matrix element when a CSF in 
        ! nI(nEl) form is provided and returns hElement of type hElement_t
        integer, intent(in) :: nI(nEl)
        HElement_t(dp) :: hel_ret
        character(*), parameter :: this_routine = "calcDiagMatEleGUGA_nI"
        
        ! have to loop over the number of spatial orbitals i , and within
        ! loop again over orbitals j > i, s indicates spatial orbitals
        integer :: iOrb, jOrb, ind, inc1, inc2, sOrb, pOrb
        real(dp) :: exchange, nOcc1, nOcc2

!         ASSERT(isProperCSF_nI(nI))

        hel_ret = ECore

        iOrb = 1
        ! loop over nI spin orbital entries: good thing is unoccupied orbitals 
        ! do not  contribute to the single matrix element part.
        do while (iOrb .le. nEl)
            ! spatial orbital index needed for get_umat_el access
            sOrb = (nI(iOrb) + 1)/2 
           ! have to check if orbital is singly or doubly occupied.
           if (isDouble(nI,iOrb)) then ! double has two part. int. contribution
               nOcc1 = 2.0_dp
               hel_ret = hel_ret + nOcc1 * GetTMatEl(nI(iOrb), nI(iOrb)) + &
                   get_umat_el(sOrb, sOrb, sOrb, sOrb)
               
               ! correctly count through spin orbitals if its a double occ.
               inc1 = 2

           else ! single occupation
               nOcc1 = 1.0_dp
               hel_ret = hel_ret + nOcc1 * GetTMatEl(nI(iOrb), nI(iOrb))
               inc1 = 1

           end if

           ! second loop:
           jOrb = iOrb + inc1
           do while (jOrb .le. nEl)
               pOrb = (nI(jOrb) + 1)/2
               ! again check for double occupancies 
               if (isDouble(nI,jOrb)) then
                   nOcc2 = 2.0_dp
                   inc2 = 2

               else 
                   nOcc2 = 1.0_dp
                   inc2 = 1

               end if
               ! standard two particle contribution
               hel_ret = hel_ret + nOcc1 * nOcc2 *( &
                   get_umat_el(sOrb,pOrb,sOrb,pOrb) - &
                   get_umat_el(sOrb,pOrb,pOrb,sOrb)/2.0_dp)

               ! calculate exchange integral part, involving Shavitt 
               ! rules for matrix elements, only contributes if both
               ! stepvectors are 1 or 2, still to be implemented..
               if (nOcc1 == 1.0_dp .and. nOcc2 == 1.0_dp) then
                   hel_ret = hel_ret - get_umat_el(sOrb,pOrb,pOrb,sOrb) * &
                       calcDiagExchangeGUGA_nI(iOrb, jOrb, nI)/2.0_dp
               end if
        
               ! increment the counters 
               jOrb = jOrb + inc2

           end do
           iOrb = iOrb + inc1
       end do

    end function calcDiagMatEleGUGA_nI

    function calcDiagMatEleGuga_ilut(ilut) result(hElement)
        ! function to calculate the diagonal matrix element if a stepvector 
        ! in ilut format is given
        integer(n_int), intent(in) :: ilut(0:niftot)
        HElement_t(dp) :: hElement

        integer :: nI(nEl)

        call decode_bit_det(nI, ilut)

        hElement = calcDiagMatEleGUGA_nI(nI)

        ! todo
    end function calcDiagMatEleGuga_ilut

    function calcDiagExchangeGUGA_nI(iOrb, jOrb, nI) result(exchange)
        ! calculates the exchange contribution to the diagonal matrix elements
        ! this is the implementation if only nI is provided
        integer, intent(in) :: iOrb, jOrb, nI(nEl)
        real(dp) :: exchange
        character(*), parameter :: this_routine = "calcDiagExchangeGUGA_nI"

        real(dp) :: bVector(nEl)
        integer :: i

        ! the b-vector is also needed for these calculations:
        bVector = calcB_vector_nI(nI)
        ! probably could use current b vector.. or reference b vector even...
        ! yes definitly. no not really since this is also used for general 
        ! diagonal matrix elements not only the current determinant in the 
        ! fciqmc loop


        exchange = 1.0_dp
        ! then i need the exchange product term between orbitals s and p
        do i = iOrb + 1, jOrb - 1
            if (.not. isDouble(nI,i)) then
                if (is_beta(nI(i))) then
                    exchange = exchange * functionA(bVector(i), 2.0_dp, 0.0_dp)&
                        * functionA(bVector(i), -1.0_dp, 1.0_dp)
                
                else
                    exchange = exchange * functionA(bVector(i), 0.0_dp, 2.0_dp) * & 
                        functionA(bVector(i), 3.0_dp, 1.0_dp)

                end if
            end if        
        end do

        if (is_beta(nI(iOrb))) then
            exchange = exchange * functionA(bVector(iOrb), 2.0_dp, 0.0_dp)

            if (is_beta(nI(jOrb))) then
                exchange = exchange * functionA(bVector(jOrb), -1.0_dp, 1.0_dp)

            else
                exchange = -exchange * functionA(bVector(jOrb), 3.0_dp, 1.0_dp)
            
            end if

        else
            exchange = exchange * functionA(bVector(iOrb), 0.0_dp, 2.0_dp)

            if (is_beta(nI(jOrb))) then
                exchange = -exchange * functionA(bVector(jOrb), -1.0_dp, 1.0_dp)

            else
                exchange = exchange * functionA(bVector(jOrb), 3.0_dp, 1.0_dp)

            end if
        end if

    end function calcDiagExchangeGUGA_nI

!     function calcDiagExchangeGUGA_ilut(ilut) result (exchange)
!         ! calculates the exchange contribution to the diagonal matrix elements
!         ! if a CSF is provided in ilut format
!         integer(n_int), intent(in) :: ilut
!         real(dp) :: exchange
! 
!         ! todo
!     end function calcDiagExchangeGUGA_ilut

    function functionA(bValue, x, y) result(r)
        ! calculated the "A" function used for Shavitts matrix element calc.
        real(dp), intent(in) :: bValue, x, y
        real(dp) :: r
        character(*), parameter :: this_routine = "functionA"

        ASSERT( bValue + y > 0.0_dp)
        ASSERT( bValue + x >= 0.0_dp)

        r = sqrt((bValue + x)/(bValue + y))

    end function functionA

end module guga_matrixElements
#endif 
