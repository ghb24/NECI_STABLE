module guga_types
    use guga_procedure_pointers, only: general_weight_dummy, general_weight_zero
    use guga_data, only: weight_data

    implicit none

    ! define a general weight obj type, to access all calculating funcitons
    type :: weight_proc
        procedure(general_weight_dummy), pointer, nopass :: minus => null()
        procedure(general_weight_dummy), pointer, nopass :: plus => null()
        procedure(general_weight_zero), pointer, nopass :: zero => null()
    end type weight_proc

    ! and store it in a general type
    type :: weight_obj
        logical :: initialized
        type(weight_data) :: dat
        type(weight_proc) :: proc
        ! is it possible to put a type into the same type? 
        ! i still only need one pointer variable, since the other one in case 
        ! of a full double excitation will be contained in the pointer weight 
        ! object as another pointer
        type(weight_obj), pointer :: ptr => null()
    end type  weight_obj

end module
