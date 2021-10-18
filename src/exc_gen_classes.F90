module exc_gen_classes
    use constants, only: dp, n_int, maxExcit, stdout
    use procedure_pointers, only: generate_excitation, gen_all_excits
    use excitation_generators, only: ExcitationGenerator_t
    use FciMCData, only: excit_gen_store_type
    use procedure_pointers, only: generate_excitation_t, generate_all_excits_t
    use bit_rep_data, only: NIfTot
    use SystemData, only: nel, tGAS
    use Determinants, only: DefDet

    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, sum
    use gasci, only: GAS_exc_gen, GAS_specification, possible_GAS_exc_gen, get_name
    use gasci_discarding, only: GAS_DiscardingGenerator_t
    use gasci_pchb, only: GAS_PCHB_ExcGenerator_t, use_supergroup_lookup, GAS_PCHB_singles_generator
    use gasci_general, only: GAS_heat_bath_ExcGenerator_t
    use gasci_disconnected, only: GAS_disc_ExcGenerator_t
    use gasci_util, only: write_GAS_info

    implicit none
    private
    public :: init_exc_gen_class, finalize_exz_gen_class, current_exc_generator, class_managed

    class(ExcitationGenerator_t), allocatable :: current_exc_generator

contains

    !> @brief
    !> This is a helper function to allow backwards compatibility.
    subroutine class_gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        call current_exc_generator%gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tParity, pGen, hel, store, part_type)
    end subroutine

    !> @brief
    !> This is a helper function to allow backwards compatibility.
    subroutine class_gen_all_excits(nI, n_excits, det_list)
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call current_exc_generator%gen_all_excits(nI, n_excits, det_list)
    end subroutine

    subroutine init_exc_gen_class()
        use SystemData, only: t_pchb_excitgen
        use pchb_excitgen, only: PCHB_FCI_excit_generator_t


        block
            if (tGAS) then
                if (GAS_exc_gen == possible_GAS_exc_gen%DISCARDING) then
                    allocate(GAS_DiscardingGenerator_t :: current_exc_generator)
                    select type(current_exc_generator)
                    type is (GAS_DiscardingGenerator_t)
                        call current_exc_generator%init(GAS_specification)
                    end select
                else if (GAS_exc_gen == possible_GAS_exc_gen%GENERAL_PCHB) then
                    allocate(GAS_PCHB_ExcGenerator_t :: current_exc_generator)
                    select type(current_exc_generator)
                    type is (GAS_PCHB_ExcGenerator_t)
                        call current_exc_generator%init(GAS_specification, use_supergroup_lookup, &
                                                        use_supergroup_lookup, GAS_PCHB_singles_generator)
                    end select
                else if (GAS_exc_gen == possible_GAS_exc_gen%GENERAL) then
                    current_exc_generator = GAS_heat_bath_ExcGenerator_t(GAS_specification)
                else if (GAS_exc_gen == possible_GAS_exc_gen%disconnected) then
                    current_exc_generator = GAS_disc_ExcGenerator_t(GAS_specification)
                end if
                write(stdout, *)
                write(stdout, '(A" is activated")') get_name(GAS_exc_gen)
                write(stdout, '(A)') 'The following GAS specification was used: '
                block
                    type(SpinProj_t) :: S_z
                    S_z = sum(calc_spin_raw(DefDet))
                    call write_GAS_info(GAS_specification, nEl, S_z, stdout)
                end block
                write(stdout, *)
                end if
        end block

        block
            if (t_pchb_excitgen) then
                allocate(PCHB_FCI_excit_generator_t :: current_exc_generator)
                select type(current_exc_generator)
                type is (PCHB_FCI_excit_generator_t)
                    call current_exc_generator%init()
                end select
            end if
        end block
    end subroutine

    subroutine class_managed(generate_excitation, gen_all_excits)
        procedure(generate_excitation_t), pointer, intent(out) :: generate_excitation
        procedure(generate_all_excits_t), pointer, intent(out) :: gen_all_excits
        generate_excitation => class_gen_exc
        gen_all_excits => class_gen_all_excits
    end subroutine

    subroutine finalize_exz_gen_class()
        if (allocated(current_exc_generator)) then
            call current_exc_generator%finalize()
            deallocate(current_exc_generator)
        end if
    end subroutine

end module
