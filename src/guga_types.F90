module guga_types
    use guga_procedure_pointers, only: general_weight_dummy, general_weight_zero
    use guga_data, only: WeightData_t

    implicit none

    private
    public :: WeightObj_t

    ! define a general weight obj type, to access all calculating funcitons
    type :: WeightProc_t
        procedure(general_weight_dummy), pointer, nopass :: minus => null()
        procedure(general_weight_dummy), pointer, nopass :: plus => null()
        procedure(general_weight_zero), pointer, nopass :: zero => null()
    end type WeightProc_t

    ! and store it in a general type
    type :: WeightObj_t
        logical :: initialized
        type(WeightData_t) :: dat
        type(WeightProc_t) :: proc
        ! is it possible to put a type into the same type?
        ! i still only need one pointer variable, since the other one in case
        ! of a full double excitation will be contained in the pointer weight
        ! object as another pointer
        type(WeightObj_t), pointer :: ptr => null()
    end type  WeightObj_t

end module
