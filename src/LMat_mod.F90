#include "macros.h"

module LMat_mod
    use constants
    use FciMCData, only: ll_node
    use HElem, only: HElement_t_SizeB
    use SystemData, only: nBasis, t12FoldSym, G1, t_mol_3_body, nel, tStoreSpinOrbs, tContact
    use MemoryManager, only: LogMemAlloc, LogMemDealloc
    use util_mod, only: get_free_unit, fuseIndex, operator(.div.)
    use gen_coul_ueg_mod, only: get_lmat_ueg, get_lmat_ua
    use shared_memory_mpi
    use sort_mod
    use hash, only: add_hash_table_entry, clear_hash_table
    use ParallelHelper, only: iProcIndex_intra
    use tc_three_body_data, only: tDampKMat, tDampLMat, tSpinCorrelator, lMatEps, &
        tSymBrokenLMat, tSparseLMat, tLMatCalc
    use procedure_pointers, only: get_lmat_el, get_lmat_el_symInternal
    use UMatCache, only: numBasisIndices  
    use LMat_aux, only: diffSpinPos, dampLMatel
    use LMat_indexing, only: lMatIndSym, lMatIndSymBroken, oldLMatInd, strideInner, strideOuter, &
        lMatIndSpin
    use LMat_calc, only: readlMatFactors, freelMatFactors, lMatCalc, lMatABCalc
    use LMat_class, only: lMat_t, sparse_lMat_t, dense_lMat_t
#ifdef USE_HDF5_
    use hdf5
#endif
    implicit none

    ! actual objects storing the 6-index integrals
    class(lMat_t), allocatable, target :: LMat, LMatAB

contains

    !------------------------------------------------------------------------------------------!
    ! Access function for six-index integrals: get a matrix element given the spin-orbitals
    !------------------------------------------------------------------------------------------!

    !------------------------------------------------------------------------------------------!

    ! this is the common logic of all 6-index lmat-acceses
    function get_lmat_el_base(a,b,c,i,j,k) result(matel)
        ! Input: a,b,c - indices of orbitals to excite to
        !        i,j,k - indices of orbitals to excite from
        ! Output: matel - matrix element of this excitation, including all exchange terms
        use UMatCache, only: gtID
        ! Gets an entry of the 3-body tensor L:
        ! L_{abc}^{ijk} - triple excitation from abc to ijk
        implicit none
        integer, value :: a,b,c
        integer, value :: i,j,k
        HElement_t(dp) :: matel
        integer(int64) :: ida, idb, idc, idi, idj, idk
        logical :: tSameSpin

        ! initialize spin-correlator check: if all spins are the same, use LMat
        ! without spin-dependent correlator, always use LMat

        ! for matrix elements involving different spins, there are three cases:
        ! each of the orbitals i,j,k can have the different spin
        ! The position is important because the spinCorrelator breaks permutation symmetry
        ! w.r.t spin -> we fix the differing spin
        tSameSpin = .not. tSpinCorrelator .or. (G1(a)%MS==G1(b)%ms .and. G1(a)%MS==G1(c)%MS)

        ! convert to spatial orbs if required
        ida = gtID(a)
        idb = gtID(b)
        idc = gtID(c)
        idi = gtID(i)
        idj = gtID(j)
        idk = gtID(k)

        matel = 0
        ! only add the contribution if the spins match

        ! TODO: Use sameSpin to determine which contributions can appear - no need for
        ! any further checks (remember to tweak for sameSpin!= possible without tSpinCorrelator)
        ! here, we add up all the exchange terms
        call addMatelContribution(i,j,k,idi,idj,idk,1)
        call addMatelContribution(j,k,i,idj,idk,idi,1)
        call addMatelContribution(k,i,j,idk,idi,idj,1)
        call addMatelContribution(j,i,k,idj,idi,idk,-1)
        call addMatelContribution(i,k,j,idi,idk,idj,-1)
        call addMatelContribution(k,j,i,idk,idj,idi,-1)
        ! if a heuristic spin-projection is done, it happens here
        if(tDampKMat .and. .not. tDampLMat) call dampLMatel(a,b,c,matel)

    contains

        subroutine addMatelContribution(p,q,r,idp,idq,idr,sgn)
            ! get a single entry of the LMat array and add it to the matrix element
            implicit none
            integer(int64), value :: idp,idq,idr
            integer, value :: p,q,r
            integer, intent(in) :: sgn
            integer(int64) :: index
            integer :: spinPos
            class(lMat_t), pointer :: lMatPtr
            real(dp) :: lMatVal
            !     integer(int64) :: ai,bj,ck

            if(G1(p)%ms == G1(a)%ms .and. G1(q)%ms == G1(b)%ms .and. G1(r)%ms == G1(c)%ms) then

                if(tContact) then
                    lMatVal = get_lmat_ueg(ida,idb,idc,idp,idq,idr)
                else if(tLMatCalc)then
                    if(tSameSpin) then
                        lMatVal = lMatCalc(ida,idb,idc,idp,idq,idr)
                    else
                        spinPos = diffSpinPos(p,q,r,a,b,c)
                        lMatVal = lMatABCalc(ida,idb,idc,idp,idq,idr, spinPos)
                    end if
                else
                    ! pick the lMat object used here according to the spin-relation
                    if(tSameSpin) then
                        lMatPtr => lMat
                        ! the indexing function is contained in the lMat object
                        index = lMatPtr%indexFunc(ida,idb,idc,idp,idq,idr)
                    else
                        ! for different spins, check which one is the different one and
                        ! call the index function accordingly
                        lMatPtr => lMatAB
                        spinPos = diffSpinPos(p,q,r,a,b,c)
                        ! the different-spin LMat assumes the first electron in the
                        ! index function has the different spin (the order of the other two does
                        ! not matter)
                        select case(spinPos)
                        case(1)
                            index = lMatPtr%indexFunc(ida,idb,idc,idp,idq,idr)
                        case(2)
                            index = lMatPtr%indexFunc(idb,ida,idc,idq,idp,idr)
                        case(3)
                            index = lMatPtr%indexFunc(idc,idb,ida,idr,idq,idp)
                        end select
                    endif
                    lMatVal = real(lMatPtr%get_elem(index),dp)
                endif
                matel = matel + sgn * lMatVal
            endif

        end subroutine addMatelContribution

    end function get_lmat_el_base

    !------------------------------------------------------------------------------------------!

    function get_lmat_el_symmetrized(a,b,c,i,j,k) result(matel)
        ! post-symmetrized access when storing non-symmetrized lMat
        implicit none
        integer, value :: a,b,c
        integer, value :: i,j,k
        HElement_t(dp) :: matel

        matel = 0.0_dp

        ! get_lmat_el_base is not symmetric in this case => symmetrize w.r.
        ! to exchange of electrons
        matel = matel + get_lmat_el_symInternal(a,b,c,i,j,k)
        matel = matel + get_lmat_el_symInternal(b,c,a,j,k,i)
        matel = matel + get_lmat_el_symInternal(c,a,b,k,i,j)
        ! + spin correction (get_lmat_el_base(ap,bp,cp,ip,jp,kp) with ap etc being the
        ! spin-swapped indices)
        matel = matel / 3.0_dp

    end function get_lmat_el_symmetrized

    !------------------------------------------------------------------------------------------!

    function get_lmat_el_spinProj(a,b,c,i,j,k) result(matel)
        ! get the spin-projected matel
        implicit none
        integer, value :: a,b,c
        integer, value :: i,j,k
        HElement_t(dp) :: matel

        ! auxiliary, spin-swapped indices
        integer :: ap, bp, cp, ip, jp, kp
        ! prefactors for the different parts
        real(dp), parameter :: directFac = 0.5625_dp
        real(dp), parameter :: swapFac = -0.1875_dp
        real(dp), parameter :: permFac = 0.0625_dp

        matel = 0.0_dp
        matel = matel + directFac * get_lmat_el_base(a,b,c,i,j,k)
        call resetAux()
        call swapSpins(ap,bp,ip,jp)
        matel = matel + swapFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
        call resetAux()
        call swapSpins(ap,cp,ip,kp)
        matel = matel + swapFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
        call resetAux()
        call swapSpins(ap,cp,ip,kp)
        call swapSpins(ap,bp,ip,jp)
        matel = matel + permFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)

    contains

        subroutine resetAux()
            implicit none
            ! reset the auxiliary variables
            ap = a
            bp = b
            cp = c
            ip = i
            jp = j
            kp = k
        end subroutine resetAux

        subroutine swapSpins(src1, src2, tgt1, tgt2)
            implicit none
            integer, intent(inout) :: src1, src2, tgt1, tgt2
            integer :: mst1, mst2, mss1, mss2

            ! get the spin values
            mss1 = mod(src1,2)
            mss2 = mod(src2,2)
            mst1 = mod(tgt1,2)
            mst2 = mod(tgt2,2)

            ! swap two spins
            src1 = src1 + mss1 - mss2
            src2 = src2 + mss2 - mss1
            tgt1 = tgt1 + mst1 - mst2
            tgt2 = tgt2 + mst2 - mst1
        end subroutine swapSpins
    end function get_lmat_el_spinProj

    !------------------------------------------------------------------------------------------!
    ! Auxiliary functions for indexing and accessing the LMat
    !------------------------------------------------------------------------------------------!

    subroutine initializeLMatPtrs()
        implicit none
        integer :: nBI

        nBI = numBasisIndices(nBasis)
        ! some typical array dimensions useful in the indexing functions
        strideInner = fuseIndex(nBI,nBI)
        strideOuter = strideInner**2

        if(tSparseLMat) then
            allocate(sparse_lMat_t :: lMat)
            allocate(sparse_lMat_t :: lMatAB)
        else
            allocate(dense_lMat_t :: lMat)
            allocate(dense_lMat_t :: lMatAB)
        endif
        ! set the LMatInd function pointer
        if(t12FoldSym) then
            lMat%indexFunc => oldLMatInd
        else if(tSymBrokenLMat) then
            ! also need to set the size of the blocks
            lMat%indexFunc => lMatIndSymBroken
        else
            lMat%indexFunc => lMatIndSym
        endif
        ! set the spin-correlator index functions
        lMatAB%indexFunc => lMatIndSpin

        ! set the get_lmat_el function pointer
        if(tSymBrokenLMat) then
            get_lmat_el => get_lmat_el_symmetrized
        else
            get_lmat_el => get_lmat_el_base
        endif

        ! the internal pointer of get_lmat_el_symmetrized
        if(tDampLMat) then
            get_lmat_el_symInternal => get_lmat_el_spinProj
        else
            get_lmat_el_symInternal => get_lmat_el_base
        endif

    end subroutine initializeLMatPtrs

    !------------------------------------------------------------------------------------------!
    ! Six-index integral I/O functions
    !------------------------------------------------------------------------------------------!

    subroutine readLMat()
        use SystemData, only: nel
        implicit none

        ! we need at least three electrons to make use of the six-index integrals
        ! => for less electrons, this can be skipped
        if(nel<=2) return

        call initializeLMatPtrs()

        if(tLMatCalc) then
            call readLMatFactors()
        else
            ! now, read lmat from file
            call lMat%read("TCDUMP","tcdump.h5")
            ! for spin-dependent LMat, also read the opp. spin matrices
            if(tSpinCorrelator) then
                ! they break permutational symmetry w.r.t spin, so the cheapest solution is
                ! to have three instances, one for each spin-permutation
                ! (the alternatives are: spin-orbitals with full symmetry or spatial orbitals
                ! without spin symmetry - both more expensive)
                call LMatAB%read("TCDUMPAB","tcdumpab.h5")
            end if
        end if
    end subroutine readLMat

    !------------------------------------------------------------------------------------------!  

    subroutine freeLMat()
        implicit none
        character(*), parameter :: t_r = "freeLMat"

        if(tLMatCalc) then
            call freeLMatFactors()
        else
            ! These are always safe to call, regardless of allocation
            call LMatAB%dealloc()
            call LMat%dealloc()
        end if
    end subroutine freeLMat

    !------------------------------------------------------------------------------------------------
    !functions for contact interaction

    function get_lmat_el_ua(a,b,c,i,j,k) result(matel)
        use SystemData, only: G1
        use UMatCache, only: gtID
        ! Gets an entry of the 3-body tensor L:
        ! L_{abc}^{ijk} - triple excitation from abc to ijk
        implicit none
        integer, value :: a,b,c
        integer :: a2,b2,c2
        integer, intent(in) :: i,j,k
        HElement_t(dp) :: matel

        ! convert to spatial orbs if required

        matel = 0

        if(G1(a)%ms == G1(b)%ms .and. G1(a)%ms.ne.G1(c)%ms ) then
            a2=a
            b2=b
            c2=c
        elseif(G1(a)%ms == G1(c)%ms .and. G1(a)%ms.ne.G1(b)%ms ) then
            a2=c
            b2=a
            c2=b
        elseif(G1(b)%ms == G1(c)%ms .and. G1(a)%ms.ne.G1(b)%ms ) then
            a2=b
            b2=c
            c2=a
        else
            return
        endif

        ! only add the contribution if the spins match
        call addMatelContribution_ua(i,j,k,1)
        call addMatelContribution_ua(j,k,i,1)
        call addMatelContribution_ua(k,i,j,1)
        call addMatelContribution_ua(j,i,k,-1)
        call addMatelContribution_ua(i,k,j,-1)
        call addMatelContribution_ua(k,j,i,-1)
    contains

        subroutine addMatelContribution_ua(p,q,r,sgn)
            implicit none
            integer, value :: p,q,r
            integer, intent(in) :: sgn
            !     integer(int64) :: ai,bj,ck

            if(G1(p)%ms == G1(a2)%ms .and. G1(q)%ms == G1(b2)%ms .and. G1(r)%ms ==G1(c2)%ms) then
                matel = matel + 2.d0 * sgn * get_lmat_ua(a2,b2,c2,p,q,r)
            endif

        end subroutine addMatelContribution_ua

    end function get_lmat_el_ua

end module LMat_mod

