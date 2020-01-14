#include "macros.h"

module sym_general_mod

    use SystemData, only: tFixLz, tNoSymGenRandExcits, iMaxLz, G1, nel, &
                          Symmetry, tKpntSym, tReltvy, t_new_real_space_hubbard, & 
                          t_tJ_model, t_heisenberg_model, nbasis, t_k_space_hubbard, &
                          t_trans_corr_hop

    use SymExcitDataMod

    use Symdata, only: nSymLabels

    use sym_mod, only: SYMPROD, RandExcitSymLabelProd, symeq

    use constants

    use DetBitOps, only: Encodebitdet, count_open_orbs

    use bit_rep_data, only: niftot, nifd

    implicit none

    interface ClassCountInd
        module procedure ClassCountInd_full_64
        module procedure ClassCountInd_full_32
        module procedure ClassCountInd_orb
    end interface

    interface ClassCountInv
        module procedure ClassCountInv_32
        module procedure ClassCountInv_64
    end interface

    interface CCIndS
        module procedure CCIndS_32
        module procedure CCIndS_64
    end interface

contains


    elemental function ClassCountInd_full_32(Spin, Sym, Mom) result(ind)

        ! Return the index into the ClassCount arrays such that variable
        ! symmetries can be easily accomodated.
        !
        ! For spin, alpha=1, beta=2; Sym = 0:nSymLabels-1; Mom = -Lmax:LMax
        ! For molecular systems, the sym is actually the symmetry of the irrep
        ! For k-points, the sym is the k-point label from SymClasses(state)
        !
        ! INTERFACED as ClassCountInd

        integer, intent(in) :: Spin, Mom
        integer(kind=int32), intent(in) :: Sym
        integer :: ind

        if(tFixLz) then
            ind = 2 * nSymLabels * (Mom + iMaxLz) + (2 * Sym + Spin)
        else
            ind = 2 * Sym + Spin
        endif

        if(tNoSymGenRandExcits) then
            if(Spin == 1) then
                ind = 1
            else
                ind = 2
            endif
        endif

    end function

    elemental function ClassCountInd_full_64(Spin, Sym, Mom) result(ind)

        ! Return the index into the ClassCount arrays such that variable
        ! symmetries can be easily accomodated.
        !
        ! For spin, alpha=1, beta=2; Sym = 0:nSymLabels-1; Mom = -Lmax:LMax
        ! For molecular systems, the sym is actually the symmetry of the irrep
        ! For k-points, the sym is the k-point label from SymClasses(state)
        !
        ! INTERFACED as ClassCountInd

        integer, intent(in) :: Spin, Mom
        integer(kind=int64), intent(in) :: Sym
        integer :: ind

        if(tFixLz) then
            ind = int(2 * nSymLabels * (Mom + iMaxLz) + (2 * Sym + Spin),sizeof_int)
        else
            ind = int(2 * Sym + Spin,sizeof_int)
        endif

        if(tNoSymGenRandExcits) then
            if(Spin == 1) then
                ind = 1
            else
                ind = 2
            endif
        endif

    end function

    elemental function ClassCountInd_orb (orb) result(ind)

        ! The same as ClassCountInd_full, only the values required are
        ! obtained for the spin orbital orb.
        !
        ! INTERFACED as ClassCountInd

        integer, intent(in) :: orb
        integer :: ind, spin, sym, mom

        ! Extract the required values
        if (is_alpha(orb)) then
            spin = 1
        else
            spin = 2
        endif

        ! This is a HACK to work around a bug in Cray Fortran v8.1.2
        if (spin == 2) spin = 2

        sym = SpinOrbSymLabel(orb)
        mom = G1(orb)%Ml

        ! To avoid cray compiler bug!
        if (spin == 2) spin = 2

        ! Calculate index as usual
        ind = ClassCountInd (spin, sym, mom)

    end function

    ! ClassCountIndex for the spatial arrays
    pure function CCIndS_32 (sym, mom) result(ind)
        integer(kind=int32), intent(in) :: sym
        integer, intent(in) :: mom
        integer :: ind

        ind =  ((ClassCountInd(1,sym,mom)-1)/2) + 1
    end function

    ! ClassCountIndex for the spatial arrays
    pure function CCIndS_64 (sym, mom) result(ind)
        integer(kind=int64), intent(in) :: sym
        integer, intent(in) :: mom
        integer :: ind

        ind =  ((ClassCountInd(1,sym,mom)-1)/2) + 1
    end function

    elemental function class_count_spin (cc_ind) result(spn)

        ! Given a class count index, return the spin of the relevant orbitals.
        ! alpha = 1, beta = 2

        integer, intent(in) :: cc_ind
        integer :: spn

        spn = 2 - mod(cc_ind, 2)

    end function

    elemental function class_count_ms (cc_ind) result(ms)

        ! Given a class count index, return 2*ms for the relevant orbiatls.

        integer, intent(in) :: cc_ind
        integer :: ms

        ms = 2 * mod(cc_ind, 2) - 1

    end function

    elemental function class_count_ml (cc_ind) result(ml)

        ! Given a class count index, return ml for the relevant orbitals

        integer, intent(in) :: cc_ind
        integer :: ml
        integer :: spn, sym2

        if (tNoSymGenRandExcits .or. .not. tFixLz) then
            ml = 0
        else
            spn = 2 - mod(cc_ind, 2)
            sym2 = (mod(cc_ind-1, 2*nSymLabels) + 1 - spn) / 2
            ml = ((cc_ind - 2*sym2 - spn) / (2 * nSymLabels)) - iMaxLz
        end if

    end function

    elemental subroutine ClassCountInv_32 (ind, sym, spin, mom)

        ! Given a Class Count Index, return the symmetry, spin and momentum
        ! of the relevant orbitals

        integer, intent(in) :: ind
        integer, intent(out) :: spin, mom
        integer(int32), intent(out) :: sym

        ! The spin is determined by the even/odd status
        ! n.b. alpha == 1, beta == 2
        spin = 2 - mod(ind, 2)

        ! How we get the symmetry/momentum depends on the parameters of the
        ! calculation
        if (tNoSymGenRandExcits) then
            mom = 0
            sym = 0
        else if (tFixLz) then
            sym = int((mod(ind-1, 2*nSymLabels)+1 - spin) / 2,int32)
            mom = int(((ind - 2 * sym - spin) / (2 * nSymLabels)) - iMaxLz,int32)
        else
            sym = int((ind - spin) / 2,int32)
            mom = 0
        end if

    end subroutine

    elemental subroutine ClassCountInv_64 (ind, sym, spin, mom)

        ! Given a Class Count Index, return the symmetry, spin and momentum
        ! of the relevant orbitals

        integer, intent(in) :: ind
        integer, intent(out) :: spin, mom
        integer(int64), intent(out) :: sym

        ! The spin is determined by the even/odd status
        ! n.b. alpha == 1, beta == 2
        spin = 2 - mod(ind, 2)

        ! How we get the symmetry/momentum depends on the parameters of the
        ! calculation
        if (tNoSymGenRandExcits) then
            mom = 0
            sym = 0
        else if (tFixLz) then
            sym = (mod(ind-1, 2*nSymLabels)+1 - spin) / 2
            mom = int(((ind - 2 * sym - spin) / (2 * nSymLabels)) - iMaxLz)
        else
            sym = (ind - spin) / 2
            mom = 0
        end if

    end subroutine

    function SymAllowedExcit(nI, nJ, ic, ex) result(bValid)

        ! Provide a check that the determinant nJ is valid according to the
        ! symmetry specification of the determinant nI. This also checks that
        ! the excitation relating these two determinants (ic, ex) is valid.

        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2, ic)
        logical :: bValid

        integer :: exLevel, ml1, ml2, i
        integer :: sym_prod_i, sym_prod_j
        type(Symmetry) :: sym_prod1, sym_prod2
        integer(n_int) :: ilut(0:niftot)

        logical :: IsValidDet ! In .F file.
        integer :: iGetExcitLevel ! In .F file

        ! Default initial value
        bValid = .true.

        ! Check reported excitation level
        exLevel = iGetExcitLevel(nI, nJ, nel)
        if (exLevel /= ic) &
            bValid = .false.

        ! Check that determinant is in increasing numerical order
        if (.not. IsValidDet(nJ, nel)) &
            bValid = .false.

        ! Check that both determinants have the same overall symmetry
        sym_prod_i = 0
        sym_prod_j = 0

        if (t_k_space_hubbard) then 
            sym_prod1 = G1(nI(1))%Sym
            sym_prod2 = G1(nJ(1))%Sym
            do i = 2, nel 
                sym_prod1 = SYMPROD(G1(nI(i))%Sym, sym_prod1)
                sym_prod2 = SYMPROD(G1(nJ(i))%Sym, sym_prod2) 
            end do

            if (.not. SYMEQ(sym_prod1, sym_prod2)) bValid = .false.
        else 
            !todo: maybe i also have to exclude the k-space hubbard case here!
            if (.not. (t_new_real_space_hubbard .or.t_tJ_model .or. t_heisenberg_model)) then 
                do i = 1, nel
                    sym_prod_i = RandExcitSymLabelProd(SymInvLabel(SpinOrbSymLabel(nI(i))), sym_prod_i)
                    sym_prod_j = RandExcitSymLabelProd(SymInvLabel(SpinOrbSymLabel(nJ(i))), sym_prod_j)
                end do
            end if
        end if

        if (sym_prod_i /= sym_prod_j) &
            bValid = .false.

        if (.not. (t_new_real_space_hubbard .or.t_tJ_model .or. t_heisenberg_model)) then 
            ! Check the symmetry properties of the excitation matrix
           if (.not. tNoSymGenRandExcits .and. .not. tKPntSym) then
              bValid = bValid .and. IsSymAllowedExcitMat(ex,ic)
            end if
        end if
        ! should i do extra tests for heisenberg and tJ? i think so 
        if (t_new_real_space_hubbard) then 
            if (t_trans_corr_hop) then 
                if (.not. (ic == 1 .or. ic == 2)) bValid = .false.
            else
                if (ic /= 1) bValid = .false. 
            end if
        end if

        if (t_tJ_model) then 
            if (nel >= nbasis/2) bValid = .false.

            call Encodebitdet(nJ, ilut)
            ! check if we have doubly occupied orbitals:
            if ((nel - count_open_orbs(ilut(0:nifd)))/2 > 0) bValid = .false.
        end if 

        if (t_heisenberg_model) then 
            if (ic /= 2) bValid = .false. 
            call Encodebitdet(nJ, ilut) 
            if (count_open_orbs(ilut) /= nbasis/2) bValid = .false. 
            if ((nel - count_open_orbs(ilut(0:nifd)))/2 > 0) bValid = .false.
        end if

        ! Check that Lz angular momentum projection is preserved if necessary
        if (tFixLz) then
            ml1 = sum(G1(ex(1, 1:ic))%ml)
            ml2 = sum(G1(ex(2, 1:ic))%ml)
            if (ml1 /= ml2) bValid = .false.
        end if

      end function SymAllowedExcit

      function IsSymAllowedExcitMat(ex,ic) result(bValid)
        integer, intent(in) :: ex(2,ic), ic
        logical :: bValid

        type(symmetry) :: sym_prod1, sym_prod2
        integer :: ms1, ms2, i
        bValid = .true.
        
        ! Check the symmetry properties of the excitation matrix
        sym_prod1 = G1(ex(1,1))%Sym
        sym_prod2 = G1(ex(2,1))%Sym
        ms1 = G1(ex(1,1))%ms
        ms2 = G1(ex(2,1))%ms              
        do i = 2, ic
           sym_prod1 = SYMPROD(sym_prod1, G1(ex(1,i))%Sym)
           sym_prod2 = SYMPROD(sym_prod2, G1(ex(2,i))%Sym)
           ms1 = ms1 + G1(ex(1,i))%ms
           ms2 = ms2 + G1(ex(2,i))%ms
        end do
        if (.not. SYMEQ(sym_prod1, sym_prod2)) bValid = .false.
        if (ms1 /= ms2 .and.(.not.tReltvy)) bValid = .false.
      end function IsSymAllowedExcitMat
  
end module
