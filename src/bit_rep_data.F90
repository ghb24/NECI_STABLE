#include "macros.h"
module bit_rep_data

    use CalcData, only: tUseRealCoeffs
    use constants

    implicit none

    ! Structure of a bit representation:

    ! | 0-NIfD: Det | Sign(Re) | Sign(Im) | Flags |
    !
    ! -------
    ! (NIfD + 1) * 64-bits              Orbital rep.
    !  1         * 32-bits              Signs (Re)
    ! (1         * 32-bits if needed)   Signs (Im)
    ! (1         * 32-bits if needed)   Flags

    ! save all the bit rep indices and lenghts in one data structure for
    ! a more clear representation
    type :: BitRep_t
        ! number of integers (-1) to store the orbital occupation
        integer :: len_orb      = -1
        ! position of the first entry for walker population
        integer :: ind_pop      =  -1
        ! length necessar to store the population (for complex and mneci..)
        integer :: len_pop      = -1
        ! the index where the flag is stored
        integer :: ind_flag     = -1
        ! the total length of the bit-representation
        integer :: len_tot      = -1
        ! the length how much of the bit-rep is broadcasted
        integer :: len_bcast    = -1
        ! the index for truncated spawning events
        integer :: ind_spawn    = -1
        ! the index for hdiag in some implementation
        integer :: ind_hdiag    = -1

        ! GUGA specific entries:
        ! the index of the communicated rdm-index (containing info about
        ! the excit-level and excit-type)
        integer :: ind_rdm_ind  = -1
        ! the index of the x0 coupling coefficient contribution
        integer :: ind_rdm_x0       = -1
        ! the index of the x1 coupling coefficient contribution
        integer :: ind_rdm_x1       = -1
        ! the x0 element (in the hamiltonian matrix element calc.)
        integer :: ind_x0       = -1
        ! the x1 element -||-
        integer :: ind_x1           = -1
        ! the delta b element in the stochastic excit-gen
        integer :: ind_b            = -10

        ! RDM specific entries:
        ! the index of the rdm biasing factor
        integer :: ind_rdm_fac      = -1
        ! the flags (initiator) of the parent
        integer :: ind_parent_flag  = -1
        ! the index where the information where the spawn came from is stored
        integer :: ind_source       = -1
        ! the index of the parent ci coeff ( in the enlarged SpawnedParts!)
        integer :: ind_parent   = -1

    end type BitRep_t

    ! make global data structure for the bit-rep indices
    ! also initialize them to -1 so we can easily spot uninitialized stuff
    type(BitRep_t) :: IlutBits = BitRep_t(), &
                      IlutBitsParent = BitRep_t(), &
                      GugaBits = BitRep_t()

    integer :: nIfTot  ! Upper bound of bit representation. In form 0:NIfTot
    integer :: nIfD    ! Final byte representing spatial/spin orbitals

    integer :: nIfTotKP ! Upper bound of krylov_vecs.

    integer :: nIfGUGA ! number of integers needed for the GUGA CSFs and flags

    ! Has the RDM component of the bit representation (For bcast) been inited.
    logical :: bit_rdm_init

    ! Flags which we can store
    integer :: flag_counter

    logical :: tuseflags = .true.

    integer, parameter :: flag_deterministic = 0, &
                          flag_determ_parent = 1, &
                          flag_trial = 2, &
                          flag_connected = 3, &
                          flag_prone = 4, &
                          flag_rescale = 5, &
                          flag_deltaB_single = 6, & ! new flags added for GUGA
                          flag_deltaB_double = 7, & ! new flags added for GUGA
                          flag_deltaB_sign = 8, &   ! new flags added for GUGA
                          flag_ic0_spawn = 9, &
                          flag_death_done = 10, &
                          flag_negative_sign = 11, &
                          flag_large_matel = 12


#ifdef PROG_NUMRUNS_
    integer, parameter :: flag_initiator(lenof_sign_max) &
                            = (/ 13, 14, 15, 16, 17, 18, 19, &
                                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32/), &
                          flag_adi_checked = 33, &
                          flag_static_init(lenof_sign_max) &
                            = (/34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, &
                                45, 46, 47, 48, 49, 50, 51, 52, 53/), &
                          flag_removed = 54, &
                          num_flags = 55
#else
    integer, parameter :: flag_initiator(2) = (/ 13, 14/), &
                          flag_adi_checked = 15, &
                          flag_static_init(2) = (/16, 17/), &
                          flag_removed = 18, &
                          num_flags = 19
#endif

contains

    pure function test_flag (ilut, flg) result(bSet)

        ! Test specified flag (0 indexed) in the bit representation.
        !
        ! In:  flg  - Integer index of flag to test
        !      ilut - Bit representation of determinant
        ! Ret: bSet - returns .true. if the flag is set, false otherwise

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
        logical :: bSet

        bSet = btest(ilut(IlutBits%ind_flag), flg)

    end function test_flag

    pure subroutine extract_sign (ilut, real_sgn)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        real(dp), intent(out) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        sgn = iLut(IlutBits%ind_pop:IlutBits%ind_pop+lenof_sign-1)
        ! transfer operates elementwise
        real_sgn = transfer(sgn, real_sgn)

    end subroutine extract_sign

end module
