module procedure_pointers

    implicit none

    !
    ! Here we define all of the procedure-pointered interfaces
    ! TODO: Package information transfer up more tightly (see SPINS branch)
    !
    abstract interface

        !
        ! Generic excitaiton generator
        subroutine generate_excitation_t(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                         ex, tParity, pGen, hel, store, part_type)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use FciMCData, only: excit_gen_store_type
            use constants
            implicit none

            integer, intent(in) :: nI(nel), exFlag
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
            integer(n_int), intent(out) :: ilutJ(0:NifTot)
            real(dp), intent(out) :: pGen
            logical, intent(out) :: tParity
            HElement_t(dp), intent(out) :: hel
            type(excit_gen_store_type), intent(inout), target :: store
            integer, intent(in), optional :: part_type

        end subroutine generate_excitation_t

        subroutine generate_single_excit_t(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use FciMCData, only: excit_gen_store_type
            use constants
            implicit none

            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(out) :: nJ(nel), ex(2, maxExcit)
            integer(n_int), intent(out) :: ilutJ(0:NIfTot)
            logical, intent(out) :: tpar
            real(dp), intent(out) :: pGen
            type(excit_gen_store_type), intent(inout), target :: store

        end subroutine generate_single_excit_t

        !
        ! Generic attempt create routine
        function attempt_create_t(nI, ilutI, wSign, nJ, ilutJ, prob, HElGen, &
                                  ic, ex, tPar, exLevel, part_type, &
                                  AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac) &
            result(child)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer, intent(in) :: nI(nel), nJ(nel), part_type, ic, ex(2, ic)
            integer, intent(in) :: exLevel
            integer(n_int), intent(in) :: ilutI(0:NifTot)
            integer(n_int), intent(inout) :: ilutJ(0:NifTot)
            real(dp), intent(in) :: wSign(lenof_sign)
            logical, intent(in) :: tPar
            real(dp), intent(inout) :: prob
            real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
            ! average number of excitations per walker for this determinant
            real(dp), intent(in) :: AvExPerWalker
            real(dp), intent(out) :: RDMBiasFacCurr
            real(dp), intent(in) :: precond_fac
            HElement_t(dp), intent(inout) :: HElGen
            real(dp) :: child(lenof_sign)

        end function

        !
        ! Generic spawning HElement routine
        function get_spawn_helement_t(nI, nJ, ilutI, ilutJ, ic, ex, tParity, &
                                      HElGen) result(hel)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer, intent(in) :: nI(nel), nJ(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
            integer, intent(in) :: ic, ex(2, ic)
            logical, intent(in) :: tParity
            HElement_t(dp), intent(in) :: HElGen
            HElement_t(dp) :: hel

        end function

        ! Generic routine to encode child determinants
        subroutine encode_child_t(ilutI, ilutJ, ic, ex)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer(n_int), intent(in) :: ilutI(0:NifTot)
            integer, intent(in) :: ic
            integer, intent(in) :: ex(2, ic)
            integer(n_int), intent(inout) :: ilutJ(0:NIfTot)

        end subroutine

        !
        ! Generic particle death routine
        function attempt_die_t(nI, Kii, wSign, exLevel, DetPosition) result(ndie)

            use SystemData, only: nel
            use constants
            implicit none

            integer, intent(in) :: nI(nel), exLevel
            real(dp), intent(in) :: Kii, wSign(lenof_sign)
            integer, intent(in), optional :: DetPosition
            real(dp), dimension(lenof_sign) :: ndie

        end function

        !
        ! Generic extract_bit_rep_avsign routine
        subroutine extract_bit_rep_avsign_t(rdm_defs, ilutI, j, nI, signI, &
                                            flagsI, IterRDMStartI, AvSignI, &
                                            store)

            ! j --> Which slot in CurrentDets are we examining.

            use rdm_data, only: rdm_definitions_t
            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use FciMCData, only: excit_gen_store_type
            use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
            use constants
            implicit none

            type(rdm_definitions_t), intent(in) :: rdm_defs
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: j
            real(dp), intent(out) :: IterRDMStartI(len_iter_occ_tot), AvSignI(len_av_sgn_tot)
            integer, intent(out) :: nI(nel), FlagsI
            real(dp), intent(out) :: SignI(lenof_sign)
            type(excit_gen_store_type), intent(inout), optional :: store

        end subroutine

        !
        ! Generic fill_rdm_diag_currdet routine
        subroutine fill_rdm_diag_currdet_t(spawn, one_rdms, ilutI, nI, ExcitLevelI, av_sign, iter_occ, tCoreSpaceDet, tLC)

            use bit_rep_data, only: NIfTot
            use constants
            use rdm_data, only: rdm_spawn_t, one_rdm_t
            use SystemData, only: nel
            implicit none

            type(rdm_spawn_t), intent(inout) :: spawn
            type(one_rdm_t), intent(inout) :: one_rdms(:)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: nI(nel), ExcitLevelI
            real(dp), intent(in) :: av_sign(:), iter_occ(:)
            logical, intent(in), optional :: tCoreSpaceDet, tLC

        end subroutine

        !
        ! Generic UMAT element routine (2e integrals)
        function get_umat_el_t(i, j, k, l) result(hel)

            use constants
            implicit none

            integer, intent(in) :: i, j, k, l
            HElement_t(dp) :: hel

        end function

        ! generic lMat element routine (3e integrals)
        function get_lmat_el_t(a, b, c, i, j, k) result(hel)
            use constants
            implicit none

            integer, value :: a, b, c
            integer, value :: i, j, k
            HElement_t(dp) :: hel

        end function get_lmat_el_t

!         subroutine generate_all_excits_t(nI, n_excits, det_list)
!             use SystemData, only: nel
!             use constants, only: n_int
!             integer, intent(in) :: nI(nel)
!             integer, intent(out) :: n_excits
!             integer(n_int), intent(out), allocatable :: det_list(:,:)
!         end subroutine generate_all_excits_t

        pure function scale_function_t(hdiag) result(Si)
            use constants
            implicit none

            real(dp), intent(in) :: hdiag
            real(dp) :: Si

        end function scale_function_t

        pure function lMatInd_t(a, b, c, i, j, k) result(index)
            use constants, only: int64
            implicit none
            integer(int64), value, intent(in) :: a, b, c ! occupied orb indices
            integer(int64), value, intent(in) :: i, j, k ! unoccupied orb
            integer(int64) :: index
        end function lMatInd_t

        pure function shift_factor_function_t(pos, run, pop) result(f)
            ! Scale facotr function for adpative shift
            ! Input: pos - position of given determinant in CurrentDets
            ! Input: run - run for which the factor is needed
            ! Input: pop - population of given determinant
            ! Output: f - scaling factor for the shift
            use constants
            implicit none

            integer, intent(in) :: pos
            integer, intent(in) :: run
            real(dp), intent(in) :: pop
            real(dp) :: f

        end function shift_factor_function_t

        !> Generate all excitations for a given determinant in the ilut Format
        subroutine generate_all_excits_t(nI, n_excits, det_list)
            use SystemData, only: nel
            use constants, only: n_int
            integer, intent(in) :: nI(nel)
            integer, intent(out) :: n_excits
            integer(n_int), intent(out), allocatable :: det_list(:, :)
        end subroutine generate_all_excits_t

    end interface

    !
    ! And here are the stored procedure pointers (for use in FCIQMC)
    procedure(generate_excitation_t), pointer :: generate_excitation => null()
    procedure(generate_excitation_t), pointer :: generate_two_body_excitation => null()
    procedure(generate_single_excit_t), pointer :: generate_single_excit => null()
    procedure(attempt_create_t), pointer :: attempt_create => null()
    procedure(get_spawn_helement_t), pointer :: get_spawn_helement => null()
    procedure(get_spawn_helement_t), pointer :: get_conn_helement => null()
    procedure(encode_child_t), pointer :: encode_child => null()
    procedure(attempt_die_t), pointer :: attempt_die => null()
    procedure(extract_bit_rep_avsign_t), pointer :: extract_bit_rep_avsign => null()
    procedure(fill_rdm_diag_currdet_t), pointer :: fill_rdm_diag_currdet => null()

    !
    ! The two UMAT (2e integral) routines. The second is only used if a
    ! 'stacking' scheme is in use (i.e. caching, memoization etc.)
    procedure(get_umat_el_t), pointer :: get_umat_el => null()
    procedure(get_umat_el_t), pointer :: get_umat_el_secondary => null()

    ! the function used to scale the walkers
    procedure(scale_function_t), pointer :: scaleFunction => null()
    ! the function used to scale the shift
    procedure(shift_factor_function_t), pointer :: shiftFactorFunction => null()

    ! the procedure to generate all determinants that are connected to a given determinant
    procedure(generate_all_excits_t), pointer :: gen_all_excits => null()

    ! the function used to scale the shift
    procedure(scale_function_t), pointer :: shiftScaleFunction => null()

end module
