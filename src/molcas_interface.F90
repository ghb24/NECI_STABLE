subroutine NECImain(NECIen)

    use constants, only : dp
    use rdm_finalising, only : RDM_energy
    use NECICore_mod, only : NECICore
    implicit none
    character(64), parameter :: int_name = 'FCIDUMP', fci_inp = 'FCINP'
    real(dp), intent (out) :: NECIen

    write(6, *) "STARTING NECI from Molcas"
    call NECICore(0, .false., .false., .false., .true., int_name, fci_inp)
    ! Once we got excited states energies we will add them here to ENER array.
    NECIen = RDM_energy

end subroutine NECImain
