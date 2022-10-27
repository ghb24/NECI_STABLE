module test_pchb_uhf_excitgen
    use fruit, only: assert_true
    use constants, only: dp
    ! use gasci_pchb_uhf, only:
    use pchb_uhf_excit, only: PCHB_UHF_FCI_excit_generator_t
        !! SUT
    use util_mod, only: near_zero
    implicit none
    private
    public :: pchb_uhf_test_driver
    ! @jph
contains

    subroutine pchb_uhf_test_driver()
        type(PCHB_UHF_FCI_excit_generator_t) :: exc_generator
        ! @jph
    end subroutine pchb_uhf_test_driver

end module test_pchb_uhf_excitgen