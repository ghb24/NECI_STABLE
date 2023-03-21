#include "macros.h"
#:include "macros.fpph"
#:set countable_types = {'integer': {'int32', 'int64'}}
#:set primitive_types = {'integer': {'int32', 'int64'}, 'real': {'sp', 'dp'}}

module sets_mod
    use constants, only: int32, int64, sp, dp
    use util_mod, only: binary_search_int, stop_all
    use sort_mod, only: sort
    use growing_buffers, only: buffer_int32_1D_t, buffer_int64_1D_t
    implicit none
    private
    public :: subset, set, is_set, is_sorted, special_union_complement, disjoint, &
              operator(.cap.), operator(.complement.), operator(.U.), &
              operator(.in.), empty_int

    integer :: empty_int(0) = [integer::]

    !> Check if V is sorted.
    interface is_sorted
    #:for T, kinds in primitive_types.items()
    #:for kind in kinds
        module procedure is_sorted_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Check if A and B are disjoint.
    interface disjoint
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure disjoint_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Create a set out of A
    interface set
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure set_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Check if a given array is a set (ordered and unique elements)
    interface is_set
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure is_set_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Check if A is a subset of B.
    interface subset
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure subset_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Calculate the union A ∪ B
    interface operator(.U.)
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure union_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Calculate the intersection A ∩ B
    interface operator(.cap.)
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure intersect_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Calculate the complement A / B
    interface operator (.complement.)
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure complement_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Check if element is contained in set.
    interface operator(.in.)
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
        module procedure test_in_${T}$_${kind}$
    #:endfor
    #:endfor
    end interface

    !> Specialiced function with assumptions that speed up performance.
    !> Merge B into A and remove values that are in C.
    !> The result can be written with set notation as A ∪ B / C.
    !> Preconditions (not tested!):
    !>      1. C is a subset of A
    !>      2. A and B are disjoint
    !>      3. B and C are disjoint
    !>      4. A, B, and C are sorted.
    !> The result will be sorted.
    interface special_union_complement
        #:for T, kinds in countable_types.items()
        #:for kind in kinds
        module procedure special_union_complement_${T}$_${kind}$
        #:endfor
        #:endfor
    end interface

contains

    #:for T, kinds in primitive_types.items()
    #:for kind in kinds
    pure function is_sorted_${T}$_${kind}$ (V, ascending) result(res)
        ${T}$ (${kind}$), intent(in) :: V(:)
        logical, intent(in), optional :: ascending
        logical :: ascending_
        logical :: res

        integer :: i

        @:def_default(ascending_, ascending, .true.)

        res = .true.
        if (ascending_) then
            do i = 1, size(V) - 1
                if (V(i) > V(i + 1)) then
                    res = .false.
                    return
                end if
            end do
        else
            do i = 1, size(V) - 1
                if (V(i) < V(i + 1)) then
                    res = .false.
                    return
                end if
            end do
        end if
    end function is_sorted_${T}$_${kind}$
    #:endfor
    #:endfor

    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function set_${T}$_${kind}$ (V) result(res)
        ${T}$ (${kind}$), intent(in) :: V(:)
        integer(${kind}$), allocatable :: res(:)
        integer(${kind}$) :: sorted(size(V)), previous
        type(buffer_${kind}$_1D_t) :: buffer
        integer :: i

        sorted = V
        call sort(sorted)

        call buffer%init()
        if (size(sorted) > 0) then
            call buffer%push_back(sorted(1))
            previous = sorted(1)
        end if
        do i = 2, size(sorted)
            if (previous /= sorted(i)) then
                call buffer%push_back(sorted(i))
                previous = sorted(i)
            end if
        end do
        call buffer%dump_reset(res)
    end function set_${T}$_${kind}$
    #:endfor
    #:endfor

    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function is_set_${T}$_${kind}$ (V) result(res)
        ${T}$ (${kind}$), intent(in) :: V(:)
        logical :: res
        integer(${kind}$) :: previous
        integer :: i

        res = is_sorted(V)
        if (res) then
            if (size(V) > 0) previous = V(1)
            do i = 2, size(V)
                if (V(i) == previous) then
                    res = .false.
                    return
                else
                    previous = V(i)
                end if
            end do
        end if
    end function is_set_${T}$_${kind}$
    #:endfor
    #:endfor

    ! check if A and B are disjoint
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function disjoint_${T}$_${kind}$ (A, B) result(res)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:)
        logical :: res
        character(*), parameter :: this_routine = "disjoint_${T}$_${kind}$"

        integer :: i, j

            @:pure_ASSERT(is_sorted(A))
            @:pure_ASSERT(is_sorted(B))

        res = .true.
        i = 1; j = 1
        do while (i <= size(A) .and. j <= size(B))
            if (A(i) < B(j)) then
                i = i + 1
            else if (A(i) > B(j)) then
                j = j + 1
            else
                res = .false.
                return
            end if
        end do
    end function disjoint_${T}$_${kind}$
    #:endfor
    #:endfor

    !> Check if A is a subset of B
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function subset_${T}$_${kind}$ (A, B) result(res)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:)
        logical :: res
        character(*), parameter :: this_routine = "subset_${T}$_${kind}$"

        integer :: i, j

            @:pure_ASSERT(is_sorted(A))
            @:pure_ASSERT(is_sorted(B))

        if (size(A) == 0) then
            res = .true.
            return
        end if

        res = .false.

        if (size(A) > size(B)) return
        if (A(1) < B(1) .or. A(size(A)) > B(size(B))) return

        i = 1; j = 1
        do while (i <= size(A))
            if (j > size(B)) then
                return
            else if (i == size(A) .and. j == size(B) .and. A(i) /= B(j)) then
                return
            else if (A(i) < B(j)) then
                i = i + 1
            else if (A(i) > B(j)) then
                j = j + 1
            else if (A(i) == B(j)) then
                i = i + 1
                j = j + 1
            end if
        end do
        res = .true.
    end function subset_${T}$_${kind}$
    #:endfor
    #:endfor

    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function special_union_complement_${T}$_${kind}$ (A, B, C) result(D)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:), C(:)
        ${T}$ (${kind}$) :: D(size(A) + size(B) - size(C))
        character(*), parameter :: this_routine = 'union_complement'

        integer :: i, j, k, l

        @:pure_ASSERT(is_sorted(A))
        @:pure_ASSERT(is_sorted(B))
        @:pure_ASSERT(is_sorted(C))
        @:pure_ASSERT(disjoint(A, B))
        @:pure_ASSERT(disjoint(B, C))
        @:pure_ASSERT(subset(C, A))

        i = 1; j = 1; k = 1; l = 1
        do while (l <= size(D))
            ! Only indices from C have to be added to A
            ! We use assumption that B is a subset of A
            if (i > size(A)) then
                D(l) = B(k)
                k = k + 1
                l = l + 1
                ! Neither indices from B have to be deleted in A
                ! nor indices from C have to be added from C to A.
            else if (j > size(C) .and. k > size(B)) then
                D(l) = A(i)
                i = i + 1
                l = l + 1
                ! No more indices from B have to be deleted in A
            else if (j > size(C)) then
                if (A(i) < B(k)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else if (A(i) > B(k)) then
                    D(l) = B(k)
                    k = k + 1
                    l = l + 1
                end if
                ! No more indices have to be added from C to A
            else if (k > size(B)) then
                if (A(i) /= C(j)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else
                    i = i + 1
                    j = j + 1
                end if
                ! Normal case:
                ! Merge C sorted into A excluding values from B.
            else if (A(i) < B(k)) then
                if (A(i) /= C(j)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else
                    i = i + 1
                    j = j + 1
                end if
            else if (A(i) > B(k)) then
                if (A(i) /= C(j)) then
                    D(l) = B(k)
                    k = k + 1
                    l = l + 1
                else
                    D(l) = B(k)
                    i = i + 1
                    j = j + 1
                    k = k + 1
                    l = l + 1
                end if
            end if
        end do
    end function special_union_complement_${T}$_${kind}$
    #:endfor
    #:endfor

    !> Return A ∪ B
    !> Assume:
    !>      1. A and B are sorted.
    !> The result will be sorted.
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function union_${T}$_${kind}$ (A, B) result(D)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:)
        ${T}$ (${kind}$), allocatable :: D(:)
        character(*), parameter :: this_routine = 'union_complement'

        ${T}$ (${kind}$), allocatable :: tmp(:)
        integer :: i, j, l

        @:pure_ASSERT(is_sorted(A))
        @:pure_ASSERT(is_sorted(B))

        allocate(tmp(size(A) + size(B)))

        i = 1; j = 1; l = 1
        do while (l <= size(tmp))
            if (i > size(A) .and. j > size(B)) then
                exit
            else if (i > size(A)) then
                tmp(l) = B(j)
                j = j + 1
                l = l + 1
            else if (j > size(B)) then
                tmp(l) = A(i)
                i = i + 1
                l = l + 1
            else if (A(i) == B(j)) then
                tmp(l) = A(i)
                i = i + 1
                j = j + 1
                l = l + 1
            else if (A(i) < B(j)) then
                tmp(l) = A(i)
                i = i + 1
                l = l + 1
            else if (A(i) > B(j)) then
                tmp(l) = B(j)
                j = j + 1
                l = l + 1
            end if
        end do

        associate(final_size => l - 1)
            D = tmp(:final_size)
        end associate
        @:pure_ASSERT(is_sorted(D))
    end function union_${T}$_${kind}$
    #:endfor
    #:endfor

    !> Return A ∩ B
    !> Assume:
    !>      1. A and B are sorted.
    !> The result will be sorted.
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function intersect_${T}$_${kind}$ (A, B) result(D)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:)
        ${T}$ (${kind}$), allocatable :: D(:)
        character(*), parameter :: this_routine = 'union_complement'

        ${T}$ (${kind}$), allocatable :: tmp(:)
        integer :: i, j, l

        @:pure_ASSERT(is_sorted(A))
        @:pure_ASSERT(is_sorted(B))

        allocate(tmp(min(size(A), size(B))))

        i = 1; j = 1; l = 1
        do while (l <= size(tmp))
            if (i > size(A) .or. j > size(B)) then
                exit
            else if (A(i) == B(j)) then
                tmp(l) = A(i)
                i = i + 1
                j = j + 1
                l = l + 1
            else if (A(i) < B(j)) then
                i = i + 1
            else if (A(i) > B(j)) then
                j = j + 1
            end if
        end do

        associate(final_size => l - 1)
            D = tmp(:final_size)
        end associate
        @:pure_ASSERT(is_sorted(D))
    end function intersect_${T}$_${kind}$
    #:endfor
    #:endfor

    !> Return A / B
    !> Assume:
    !>      1. A and B are sorted.
    !> The result will be sorted.
    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function complement_${T}$_${kind}$ (A, B) result(D)
        ${T}$ (${kind}$), intent(in) :: A(:), B(:)
        ${T}$ (${kind}$), allocatable :: D(:)
        character(*), parameter :: this_routine = 'union_complement'

        ${T}$ (${kind}$), allocatable :: tmp(:)
        integer :: i, j, l

        @:pure_ASSERT(is_sorted(A))
        @:pure_ASSERT(is_sorted(B))

        allocate(tmp(size(A)))

        i = 1; j = 1; l = 1
        do while (l <= size(tmp))
            if (i > size(A)) then
                exit
            else if (j > size(B)) then
                tmp(l) = A(i)
                i = i + 1
                l = l + 1
            else if (A(i) == B(j)) then
                i = i + 1
                j = j + 1
            else if (A(i) < B(j)) then
                tmp(l) = A(i)
                i = i + 1
                l = l + 1
            else if (A(i) > B(j)) then
                j = j + 1
            end if
        end do

        associate(final_size => l - 1)
            D = tmp(:final_size)
        end associate
        @:pure_ASSERT(is_sorted(D))
    end function complement_${T}$_${kind}$
    #:endfor
    #:endfor

    #:for T, kinds in countable_types.items()
    #:for kind in kinds
    pure function test_in_${T}$_${kind}$ (element, set) result(res)
        ${T}$ (${kind}$), intent(in) :: element, set(:)
        character(*), parameter :: this_routine = 'test_in_${T}$_${kind}$'
        logical :: res
        @:pure_ASSERT(is_sorted(set))
        res = binary_search_int(set, element) /= -1
    end function
    #:endfor
    #:endfor

end module sets_mod
