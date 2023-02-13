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
    use MPI_wrapper, only: iProcIndex_intra
    use tc_three_body_data, only: lMatEps, tSparseLMat, tLMatCalc, tSymBrokenLMat, tHDF5LMat
    use UMatCache, only: numBasisIndices
    use LMat_indexing, only: lMatIndSym, lMatIndSymBroken, oldLMatInd, strideInner, strideOuter, &
                             lMatIndSpin
    use LMat_calc, only: readlMatFactors, freelMatFactors, lMatCalc, &
                         read_rs_lmat_factors
    use LMat_class, only: lMat_t, sparse_lMat_t, dense_lMat_t
#ifdef USE_HDF5_
    use hdf5
#endif
    use IntegralsData, only: t_use_tchint_lib, tchint_mode, t_rs_factors
#ifdef USE_TCHINT_
    use tchint, only: tc_matel, tchint_init
#endif
    implicit none

    ! actual objects storing the 6-index integrals
    class(lMat_t), allocatable, target :: LMat

contains

    !------------------------------------------------------------------------------------------!
    ! Access function for six-index integrals: get a matrix element given the spin-orbitals
    !------------------------------------------------------------------------------------------!

    !------------------------------------------------------------------------------------------!

    ! this is the common logic of all 6-index lmat-acceses
    function get_lmat_el(a, b, c, i, j, k) result(matel)
        ! Input: a,b,c - indices of orbitals to excite to
        !        i,j,k - indices of orbitals to excite from
        ! Output: matel - matrix element of this excitation, including all exchange terms
        use UMatCache, only: gtID
        use SystemData, only: t_exclude_pure_parallel
        ! Gets an entry of the 3-body tensor L:
        ! L_{abc}^{ijk} - triple excitation from abc to ijk
        integer, value :: a, b, c
        integer, value :: i, j, k
        HElement_t(dp) :: matel
        integer(int64) :: ida, idb, idc, idi, idj, idk

        ! initialize spin-correlator check: if all spins are the same, use LMat
        ! without spin-dependent correlator, always use LMat

        matel = h_cast(0.0_dp)

        if (t_exclude_pure_parallel) then
            if (     (G1(a)%ms == G1(b)%ms .and. G1(a)%ms == G1(c)%ms) &
                .or. (G1(i)%ms == G1(j)%ms .and. G1(i)%ms == G1(k)%ms)) then
                return
            end if
        end if

        ! convert to spatial orbs if required
        ida = gtID(a)
        idb = gtID(b)
        idc = gtID(c)
        idi = gtID(i)
        idj = gtID(j)
        idk = gtID(k)

        ! only add the contribution if the spins match

        ! here, we add up all the exchange terms
        call addMatelContribution(i, j, k, idi, idj, idk, 1)
        call addMatelContribution(j, k, i, idj, idk, idi, 1)
        call addMatelContribution(k, i, j, idk, idi, idj, 1)
        call addMatelContribution(j, i, k, idj, idi, idk, -1)
        call addMatelContribution(i, k, j, idi, idk, idj, -1)
        call addMatelContribution(k, j, i, idk, idj, idi, -1)

    contains

        subroutine addMatelContribution(p, q, r, idp, idq, idr, sgn)
            ! get a single entry of the LMat array and add it to the matrix element
            integer(int64), value :: idp, idq, idr
            integer, value :: p, q, r
            integer, intent(in) :: sgn
            integer(int64) :: index
            integer :: spinPos
            class(lMat_t), pointer :: lMatPtr
            real(dp) :: lMatVal
            !     integer(int64) :: ai,bj,ck

            if (G1(p)%ms == G1(a)%ms .and. G1(q)%ms == G1(b)%ms .and. G1(r)%ms == G1(c)%ms) then

                if (tContact) then
                    lMatVal = get_lmat_ueg(ida, idb, idc, idp, idq, idr)
                else if (tLMatCalc .or. t_rs_factors) then
                    lMatVal = lMatCalc(ida, idb, idc, idp, idq, idr)

                else
                    ! the indexing function is contained in the lMat object
                    index = lMat%indexFunc(ida, idb, idc, idp, idq, idr)
                    lMatVal = real(lMat%get_elem(index), dp)
                end if
                matel = matel + sgn * lMatVal
            end if

        end subroutine addMatelContribution

    end function get_lmat_el

    !------------------------------------------------------------------------------------------!
    ! Auxiliary functions for indexing and accessing the LMat
    !------------------------------------------------------------------------------------------!

    subroutine initializeLMatPtrs()
        integer :: nBI

        nBI = numBasisIndices(nBasis)
        ! some typical array dimensions useful in the indexing functions
        strideInner = fuseIndex(nBI, nBI)
        strideOuter = strideInner**2

        if(allocated(lMat)) deallocate(lMat)
        if (tSparseLMat) then
            allocate(sparse_lMat_t :: lMat)
        else
            allocate(dense_lMat_t :: lMat)
        end if
        ! set the LMatInd function pointer
        if (t12FoldSym) then
            lMat%indexFunc => oldLMatInd
        else if (tSymBrokenLMat) then
            ! also need to set the size of the blocks
            lMat%indexFunc => lMatIndSymBroken
        else
            lMat%indexFunc => lMatIndSym
        end if

    end subroutine initializeLMatPtrs

    !------------------------------------------------------------------------------------------!
    ! Six-index integral I/O functions
    !------------------------------------------------------------------------------------------!

    subroutine setup_tchint_ints()
#ifndef USE_TCHINT_
      character(*), parameter :: t_r = "setup_tchint_ints"
#endif

      if(t_use_tchint_lib) then
#ifdef USE_TCHINT_
        call tchint_init()
#else
        call stop_all(t_r, "Did not compile with TCHINT support")
#endif
      end if
    end subroutine setup_tchint_ints

    subroutine readLMat()
        use SystemData, only: nel
        character(255) :: tcdump_name
        ! we need at least three electrons to make use of the six-index integrals
        ! => for less electrons, this can be skipped
        if (nel <= 2) return

        if(t_use_tchint_lib) return

        call initializeLMatPtrs()

        if (tLMatCalc) then
            call readLMatFactors()
        else if (t_rs_factors) then
            call read_rs_lmat_factors()
        else
            ! now, read lmat from file
            if (tHDF5LMat) then
                tcdump_name = "tcdump.h5"
            else
                tcdump_name = "TCDUMP"
            end if
            call lMat%read(trim(tcdump_name))
        end if
    end subroutine readLMat

    !------------------------------------------------------------------------------------------!

    subroutine freeLMat()
#ifndef USE_TCHINT_
        character(*), parameter :: t_r = "freeLMat"
#endif

        if (t_use_tchint_lib) then
#ifdef USE_TCHINT_
            call tchint_finalize()
#else
            call stop_all(t_r, "Did not compile with TCHINT support")
#endif
        else
            if (tLMatCalc .or. t_rs_factors) then
                call freeLMatFactors()
            else
                call LMat%safe_dealloc()
            end if
        endif
    end subroutine freeLMat

    !------------------------------------------------------------------------------------------------!

    pure function external_lMat_matel(nI, ex) result(matel)
      integer, intent(in) :: nI(:)
      integer, intent(in) :: ex(:,:)
      HElement_t(dp) :: matel

#ifdef USE_TCHINT_
      matel = tc_matel(nI, ex)
#else
      character(*), parameter :: t_r = "external_lMat_matel"
      unused_var(nI); unused_var(ex); unused_var(matel)
      call stop_all(t_r, "Did not compile with TCHINT support")
#endif

    end function external_lMat_matel

    !------------------------------------------------------------------------------------------------
    !functions for contact interaction

    pure function get_lmat_el_ua(a, b, c, i, j, k) result(matel)
        use SystemData, only: G1
        use UMatCache, only: gtID
        ! Gets an entry of the 3-body tensor L:
        ! L_{abc}^{ijk} - triple excitation from abc to ijk
        integer, value :: a, b, c
        integer :: a2, b2, c2
        integer, intent(in) :: i, j, k
        HElement_t(dp) :: matel

        ! convert to spatial orbs if required

        matel = 0

        if (G1(a)%ms == G1(b)%ms .and. G1(a)%ms /= G1(c)%ms) then
            a2 = a
            b2 = b
            c2 = c
        else if (G1(a)%ms == G1(c)%ms .and. G1(a)%ms /= G1(b)%ms) then
            a2 = c
            b2 = a
            c2 = b
        else if (G1(b)%ms == G1(c)%ms .and. G1(a)%ms /= G1(b)%ms) then
            a2 = b
            b2 = c
            c2 = a
        else
            return
        end if

        ! only add the contribution if the spins match
        call addMatelContribution_ua(i, j, k, 1, matel)
        call addMatelContribution_ua(j, k, i, 1, matel)
        call addMatelContribution_ua(k, i, j, 1, matel)
        call addMatelContribution_ua(j, i, k, -1, matel)
        call addMatelContribution_ua(i, k, j, -1, matel)
        call addMatelContribution_ua(k, j, i, -1, matel)
    contains

        pure subroutine addMatelContribution_ua(p, q, r, sgn, matel)
            integer, value :: p, q, r
            integer, intent(in) :: sgn
            HElement_t(dp), intent(inout) :: matel

            if (G1(p)%ms == G1(a2)%ms .and. G1(q)%ms == G1(b2)%ms .and. G1(r)%ms == G1(c2)%ms) then
                matel = matel + 2.d0 * sgn * get_lmat_ua(a2, b2, c2, p, q, r)
            end if

        end subroutine addMatelContribution_ua

    end function get_lmat_el_ua

end module LMat_mod
