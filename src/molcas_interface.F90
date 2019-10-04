subroutine NECImain(fcidmp, input_name, NECIen)

    use constants, only : dp, iout
    use rdm_finalising, only : RDM_energy
    use NECICore_mod, only : NECICore
    implicit none
    character(*), intent(in) :: fcidmp, input_name
    real(dp), intent (out) :: NECIen

    write(iout, *) "STARTING NECI from Molcas"
    call NECICore(0, .false., .false., .false., .true., fcidmp, input_name)
    ! Once we got excited states energies we will add them here to ENER array.
    NECIen = RDM_energy

end subroutine NECImain
