module procedure_pointers

    implicit none

    !
    ! Here we define all of the procedure-pointered interfaces
    ! TODO: Package information transfer up more tightly (see SPINS branch)
    !
    abstract interface

        !
        ! Generic excitaiton generator
        subroutine generate_excitation_t (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                          ex, tParity, pGen, hel, store, part_type)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use FciMCData, only: excit_gen_store_type
            use constants
            implicit none

            integer, intent(in) :: nI(nel), exFlag
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(out) :: nJ(nel), ic, ex(2,2)
            integer(n_int), intent(out) :: ilutJ(0:NifTot)
            real(dp), intent(out) :: pGen
            logical, intent(out) :: tParity
            HElement_t(dp), intent(out) :: hel
            type(excit_gen_store_type), intent(inout), target :: store
            integer, intent(in), optional :: part_type

        end subroutine


        !
        ! Generic attempt create routine
        function attempt_create_t (nI, ilutI, wSign, nJ, ilutJ, prob, HElGen, &
                                   ic, ex, tPar, exLevel, part_type, &
                                   AvSignCurr, RDMBiasFacCurr, precond_fac) result(child)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer, intent(in) :: nI(nel), nJ(nel), part_type, ic, ex(2,2)
            integer, intent(in) :: exLevel
            integer(n_int), intent(in) :: ilutI(0:NifTot)
            integer(n_int), intent(inout) :: ilutJ(0:NifTot)
            real(dp), intent(in) :: wSign(lenof_sign)
            logical, intent(in) :: tPar
            real(dp), intent(inout) :: prob
            real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
            real(dp), intent(out) :: RDMBiasFacCurr
            real(dp), intent(in) :: precond_fac
            HElement_t(dp), intent(in) :: HElGen
            real(dp) :: child(lenof_sign)

        end function


        !
        ! Generic spawning HElement routine
        function get_spawn_helement_t (nI, nJ, ilutI, ilutJ, ic, ex, tParity, &
                                       HElGen) result(hel)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer, intent(in) :: nI(nel), nJ(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
            integer, intent(in) :: ic, ex(2,2)
            logical, intent(in) :: tParity
            HElement_t(dp), intent(in) :: HElGen
            HElement_t(dp) :: hel

        end function


        ! Generic routine to encode child determinants
        subroutine encode_child_t (ilutI, ilutJ, ic, ex)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            implicit none

            integer(n_int), intent(in) :: ilutI(0:NifTot)
            integer, intent(in) :: ic, ex(2,2)
            integer(n_int), intent(inout) :: ilutJ(0:NIfTot)

        end subroutine


        ! Generic routine to deal with new particle statistics
        subroutine new_child_stats_t (iter_data, ilutI, nJ, ilutJ, ic, &
                                      walkExLevel, child, parent_flags, &
                                      part_type)

            use SystemData, only: nel
            use bit_rep_data, only: NIfTot
            use constants
            use FciMCData, only: fcimc_iter_data
            implicit none

            integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
            integer, intent(in) :: part_type
            real(dp), intent(in) :: child(lenof_sign)
            integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
            type(fcimc_iter_data), intent(inout) :: iter_data

        end subroutine


        !
        ! Generic particle death routine
        function attempt_die_t (nI, Kii, wSign, exLevel) result(ndie)

            use SystemData, only: nel
            use constants
            implicit none

            integer, intent(in) :: nI(nel), exLevel
            real(dp), intent(in) :: Kii, wSign(lenof_sign)
            real(dp), dimension(lenof_sign) :: ndie

        end function


        !
        ! Generic extract_bit_rep_avsign routine
        subroutine extract_bit_rep_avsign_t (rdm_defs, ilutI, j, nI, signI, &
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
        subroutine fill_rdm_diag_currdet_old_t (rdm, one_rdm, irdm, ilutI, nI, j, ExcitLevelI, tCoreSpaceDet)

            ! j --> Which slot in CurrentDets are we examining.

            use bit_rep_data, only: NIfTot
            use constants
            use rdm_data, only: one_rdm_t
            use SystemData, only: nel
            implicit none

            type(one_rdm_t), intent(inout) :: one_rdm
            integer, intent(in) :: irdm
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: nI(nel), ExcitLevelI, j
            logical, intent(in), optional :: tCoreSpaceDet

        end subroutine

        !
        ! Generic fill_rdm_diag_currdet routine
        subroutine fill_rdm_diag_currdet_t (spawn, one_rdms, ilutI, nI, ExcitLevelI, av_sign, iter_occ, tCoreSpaceDet)

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
            logical, intent(in), optional :: tCoreSpaceDet

        end subroutine


        !
        ! Generic UMAT element routine (2e integrals)
        function get_umat_el_t (i, j, k, l) result(hel)

            use constants
            implicit none

            integer, intent(in) :: i, j, k, l
            HElement_t(dp) :: hel

        end function

        function scale_function_t(hdiag) result(Si)
          use constants
          implicit none

          real(dp), intent(in) :: hdiag
          real(dp) :: Si

        end function scale_function_t


    end interface


    !
    ! And here are the stored procedure pointers (for use in FCIQMC)
    procedure(generate_excitation_t), pointer :: generate_excitation
    procedure(attempt_create_t), pointer :: attempt_create
    procedure(get_spawn_helement_t), pointer :: get_spawn_helement
    procedure(get_spawn_helement_t), pointer :: get_conn_helement
    procedure(encode_child_t), pointer :: encode_child
    procedure(new_child_stats_t), pointer :: new_child_stats
    procedure(attempt_die_t), pointer :: attempt_die
    procedure(extract_bit_rep_avsign_t), pointer :: extract_bit_rep_avsign
    procedure(fill_rdm_diag_currdet_t), pointer :: fill_rdm_diag_currdet


    !
    ! The two UMAT (2e integral) routines. The second is only used if a
    ! 'stacking' scheme is in use (i.e. caching, memoization etc.)
    procedure(get_umat_el_t), pointer :: get_umat_el
    procedure(get_umat_el_t), pointer :: get_umat_el_secondary

    ! the function used to scale the walkers
    procedure(scale_function_t), pointer :: scaleFunction

end module
