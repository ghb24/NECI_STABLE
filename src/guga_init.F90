#include "macros.h"
! guga module containing all necessary functionality needed to initialize 
! a guga simulation
module guga_init
    ! module use statements
    use SystemData, only: tCSF, tSPN, tHPHF, lNoSymmetry, STOT, nEl, &
        nBasis, tGUGA, tNoBrillouin, tExactSizeSpace, tUHF
    use hist_data, only: tHistSpawn
    use LoggingData, only: tCalcFCIMCPsi, tPrintOrbOcc
    use spin_project, only: tSpinProject
    use bit_rep_data, only: tUseFlags
    use guga_data, only: init_guga
    ! variable declaration
    implicit none

contains

    subroutine SysInitGUGA()
        ! in general this just prints out information of the state of the 
        ! guga calculation and to avoid having repeating code elements, does 
        ! not check input options and values again
     
        ! have to ensure to use flags to be able to encode deltaB values
        ! in the ilut bit representation for the excitation generation

        ! initialize the pointers...
        call init_guga()

    end subroutine SysInitGUGA
end module guga_init





         
