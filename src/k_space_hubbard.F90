#include "macros.h" 

! finally write a specific routine for the k-space hubbard, where every 
! necessary functionality for the k-space/momentum space hubbard is brought 
! together. since i want to implement a transcorrelation factor for the 
! k-space hubbard, which would necessitate a cumulative sum exciation creation 
! i decided it is better to start a new branch, instead of hacking into the 
! old branches of the code. 
! in the end i also want to combine this fully with the lattice_mod implementation
! so i can easily decide which lattice to choose and then decide if we want to 
! use a real or momentum space basis (and in the future maybe even wavelets) 
module k_space_hubbard 
    use SystemData, only: t_lattice_model, t_k_space_hubbard, t_trans_corr, & 
                    trans_corr_param, t_trans_corr_2body, trans_corr_param_2body, & 
                    nel
    use lattice_mod, only: get_helement_lattice_ex_mat, get_helement_lattice_general
    use procedure_pointers, only: get_umat_el
    use gen_coul_ueg_mod, only: get_hub_umat_el
    use constants, only: n_int, dp

    implicit none 

    ! i especially need an interface for the matrix element calculation to 
    ! implement the transcorrelated hamiltonian 
    interface get_helement_k_space_hub
        module procedure get_helement_k_space_hub_ex_mat
        module procedure get_helement_k_space_hub_general
    end interface get_helement_k_space_hub

contains 

    subroutine init_k_space_hubbard() 

        print *, " new k-space hubbard implementation init:" 

        call check_k_space_hubbard_input()

        get_umat_el => get_hub_umat_el

        call init_get_helement_k_space_hub()

    end subroutine init_k_space_hubbard

    subroutine check_k_space_hubbard_input()

        print *, "checking input for k-space hubbard:" 

        print *, "input is fine!"

    end subroutine check_k_space_hubbard_input

    subroutine init_get_helement_k_space_hub
        get_helement_lattice_ex_mat => get_helement_k_space_hub_ex_mat
        get_helement_lattice_general => get_helement_k_space_hub_general
        ! maybe i have to initialize more here, especially if we are using the 
        ! HPHF keyword I guess.. 

    end subroutine init_get_helement_k_space_hub

    function get_helement_k_space_hub_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ic, ex(2,2)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel
#ifdef __DEBUG 
        character(*), parameter :: this_routine ="get_helement_k_space_hub_ex_mat"
#endif

    end function get_helement_k_space_hub_ex_mat

    function get_helement_k_space_hub_general(nI, nJ, ic_ret) result(hel) 
        integer, intent(in) :: nI(nel), nJ(nel) 
        integer, intent(inout), optional :: ic_ret
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_helement_k_space_hub_general"
#endif

    end function get_helement_k_space_hub_general

end module k_space_hubbard
