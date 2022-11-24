module test_pchb_uhf_excitgen
    use fruit, only: assert_true
    use constants, only: dp
    ! use gasci_pchb_uhf, only:
    use pchb_uhf_excit, only: PCHB_UHF_FCI_excit_generator_t
        !! SUT
    use util_mod, only: near_zero
    use orb_idx_mod, only: calc_spin_raw
    use unit_test_helper_excitgen, only: generate_random_integrals, &
                                    init_excitgen_test, finalize_excitgen_test, &
                                    FciDumpWriter_t, test_excitation_generator
    ! TODO: I suspect (but am not confident) the imports here will be removed
    use gasci_pc_select_particles, only: PCHB_particle_selections, possible_PCHB_singles
    implicit none
    private
    public :: pchb_uhf_test_driver
    ! @jph
contains
    ! @jph need:
    ! - [x] PCHB_particle_selections
    ! - [x] PCHB_FCI_excit_generator_t
    ! - [x] possible_PCHB_singles

    subroutine pchb_uhf_test_driver()
        integer, parameter :: n_iters = 5 * 10**7
        type(PCHB_UHF_FCI_excit_generator_t) :: exc_generator
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10], n_spat_orb = 10
        logical :: successful
        ! @jph
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))

        ! TODO this needs different particle selections based on UHF I think
        call exc_generator%init(PCHB_particle_selections%UNIFORM, possible_PCHB_singles%UNIFORM)


    contains

        subroutine random_fcidump(iunit, hermitian)
            integer, intent(in) :: iunit
            logical, intent(in) :: hermitian

            call generate_random_integrals(iunit, n_el=size(det_I), n_spat_orb, &
                    sparse=0.7_dp, sparseT=0.7_dp, total_ms=calc_spin_raw(det_I), &
                    uhf=.true., hermitian=hermitian)

        end subroutine random_fcidump
    end subroutine pchb_uhf_test_driver

end module test_pchb_uhf_excitgen