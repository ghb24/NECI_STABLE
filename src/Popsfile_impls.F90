#include "macros.h"
submodule (PopsfileMod) Popsfile_impls
    use fcimc_initialisation, only: FDet
    use fcimcdata, only: fcimcstats_unit2, initiatorstats_unit, &
        complexstats_unit, EXLEVELStats_unit, fcimcstats_unit, &
        tDebug
    use LoggingData, only: tLogComplexPops, tLogEXLEVELStats
    use fcimc_initialisation, only: DeallocFCIMCMemPar, InitFCIMCCalcPar, &
        SetupParameters
    better_implicit_none

contains

    ! This routine will change the reference determinant to DetCurr. It will
    ! also re-zero all the energy estimators, since they now correspond to
    ! projection onto a different determinant.
    !
    ! n.b. NOT MODULARISED. This is a little evil, but there is an unbreakable
    !      circular dependency otherwise.
    !
    ! **** See interface in Popsfile.F90 ****
    module subroutine ChangeRefDet(DetCurr)
        INTEGER, intent(in) :: DetCurr(NEl)
        integer :: i

        FDet(1 : nEL) = DetCurr(1 : nEl)

        write(stdout, "(A)") "*** Changing the reference determinant ***"
        write(stdout, "(A)") "Switching reference and zeroing energy counters - restarting simulation"
    !
    !Initialise variables for calculation on each node
        Iter = 1
        CALL DeallocFCIMCMemPar()
        IF (iProcIndex == Root) THEN
            close(fcimcstats_unit)
            if (inum_runs == 2) close(fcimcstats_unit2)
            IF (tTruncInitiator) close(initiatorstats_unit)
            IF (tLogComplexPops) close(complexstats_unit)
            if (tLogEXLEVELStats) close(EXLEVELStats_unit)
        end if
        IF (TDebug) close(11)
        CALL SetupParameters()
        CALL InitFCIMCCalcPar()

    end subroutine ChangeRefDet

end submodule
