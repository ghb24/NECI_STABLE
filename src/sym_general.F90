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

    end functioN

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

    end functioN

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
        sym = SpinOrbSymLabel(orb)
        mom = G1(orb)%Ml

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


end module
