! The data in this module only relates to exact and Lanczos-based calculations
! of spectra, not to KP-FCIQMC calculations.

module spectral_data

    use constants, only: dp
    use FciMCData, only: hamiltonian, perturbation

    real(dp), allocatable :: pert_ground_left(:)
    real(dp), allocatable :: pert_ground_right(:)

    type(perturbation), allocatable :: left_perturb_spectral(:)
    type(perturbation), allocatable :: right_perturb_spectral(:)

    real(dp) :: left_pert_norm
    real(dp) :: right_pert_norm

    real(dp), allocatable :: trans_amps_left(:)
    real(dp), allocatable :: trans_amps_right(:)

    integer :: nomega_spectral
    real(dp) :: spectral_broadening
    real(dp) :: delta_omega_spectral
    real(dp) :: min_omega_spectral

    logical :: tIncludeGroundSpectral
    real(dp) :: spectral_ground_energy

end module spectral_data
