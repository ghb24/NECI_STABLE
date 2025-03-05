#include "macros.h"

module block_spawns

    use SystemData, only: nel
    use constants, only: n_int, stdout
    use util_mod, only: stop_all, get_free_unit
    use input_parser_mod, only: TokenIterator_t, FileReader_t, parse_definedet,&
                                ManagingFileReader_t
    use DetBitOps, only: EncodeBitDet
    use bit_rep_data, only: nifd, niftot
    use bit_reps, only: decode_bit_det
    use FciMCData, only: ll_node
    use hash, only: init_hash_table, fill_in_hash_table, clear_hash_table, &
                    hash_table_lookup
    use growing_buffers, only: buffer_int64_2D_t

    better_implicit_none

    private

    public :: init_block_list, tBlockSpawns, blocklist_name, block_list, &
              is_excitation_aborted, blocklist_hash_table

    logical :: tBlockSpawns
    character(len=255) :: blocklist_name
    integer(kind=n_int), allocatable :: block_list(:, :)
    type(ll_node), pointer :: blocklist_hash_table(:)

contains


    subroutine init_block_list(blocklistname, blocklist, blocklist_ht)
        integer(kind=n_int), allocatable, intent(out) :: blocklist(:,:)
        type(ll_node), pointer, intent(out) :: blocklist_ht(:)
        character(len=255), intent(in) :: blocklistname
        integer :: nI_blocked(nel)
        integer(kind=n_int) :: ilut_blocked(0:niftot)
        type(TokenIterator_t) :: tokens
        class(FileReader_t), allocatable :: filereader
        type(buffer_int64_2D_t) :: buffer

        call buffer%init(size(ilut_blocked))

        filereader = ManagingFileReader_t(blocklistname)

        do while (filereader%nextline(tokens, skip_empty=.false.))
            call parse_definedet(tokens, nI_blocked)

            call EncodeBitDet(nI_blocked, ilut_blocked)

            call buffer%push_back(ilut_blocked)

#ifdef DEBUG_
            write(stdout, *) "Determinant/CSF to be blocked:", nI_blocked
#endif
        end do

        call filereader%close()

        call buffer%dump_reset(blocklist)

        write(stdout, *) 'Read in', size(blocklist, 2), 'determinants/CSFs to be blocked.'

        call init_blocklist_hashtable(blocklist, blocklist_ht)

    end subroutine init_block_list

    subroutine init_blocklist_hashtable(blocklist, blocklist_ht)
        integer(kind=n_int), intent(in) :: blocklist(0:,:)
        type(ll_node), pointer, intent(out) :: blocklist_ht(:)

        allocate(blocklist_ht(size(blocklist, 2)))

        call init_hash_table(blocklist_ht)

        call fill_in_hash_table(&
            blocklist_ht, size(blocklist, 2), blocklist, size(blocklist, 2), &
            .false.)

    end subroutine init_blocklist_hashtable

    function is_excitation_aborted(excited_nJ, excited_ilut) result(tAbort)
        integer, intent(in) :: excited_nJ(nel)
        integer(kind=n_int), intent(in) :: excited_ilut(0:niftot)
        logical :: tAbort
        integer :: idx, DetHash

        tAbort = .false.

        ! At least, `gen_excit_rs_hubbard` sets nJ(1) = 0 if the excitation is 
        ! invalied
        if (excited_nJ(1) == 0) return

        if (allocated(block_list)) then
            call hash_table_lookup(&
                excited_nJ, excited_ilut, nifd, blocklist_hash_table, &
                block_list, idx, DetHash, tAbort)
        end if

    end function is_excitation_aborted

end module block_spawns
