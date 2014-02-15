#include "macros.h"

module sym_general_mod

    use SystemData, only: tFixLz, tNoSymGenRandExcits, iMaxLz, G1
    use SymExcitDataMod
    use Symdata, only: nSymLabels
    use constants

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
            sym = (mod(ind-1, 2*nSymLabels)+1 - spin) / 2
            mom = ((ind - 2 * sym - spin) / (2 * nSymLabels)) - iMaxLz
        else
            sym = (ind - spin) / 2
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
            mom = ((ind - 2 * sym - spin) / (2 * nSymLabels)) - iMaxLz
        else
            sym = (ind - spin) / 2
            mom = 0
        end if

    end subroutine


end module
