module quicksort
    use util_mod, only: arr_gt, arr_lt, swap
    use mt95, only: genrand_real2
    implicit none

contains


    subroutine test_sort_arr_int ()
        integer :: arr(3,40), i,j
        real*8 :: r

        do i=1,size(arr(:,1))
            do j=1,size(arr(1,:))
                call genrand_real2(r)
                arr(i,j) = int(r * 100)
            enddo
        enddo

        print*, 'INIT'
        do i=1,size(arr(1,:))
            print*, arr(:,i)
        enddo

        call sort_arr_int (arr)

        print*, 'FINAL'
        do i=1,size(arr(1,:))
            print*, arr(:,i)
        enddo
        call flush(6)

        call stop_all ("end", "test")
    end subroutine

    subroutine sort_arr_int (arr)

        !/* Perform a quicksort on an array of integers, arr(n). Uses the 
        !   sample code in NumericalRecipies as a base. */

        integer, intent(inout) :: arr(:,:)
        integer, parameter :: nInsertionSort = 7
        !/* sort needs auxiliary storage of length 2*log_2(n). */
        integer, parameter :: nStackMax = 50

        integer :: pivot, lo, hi, n, i, j
        !/* Oh, how lovely it would be to be able to use push/pop and not need
        !   to guess a size of the stack to start with */
        integer :: stack(nStackMax), nstack
        integer :: a(size(arr(:,1)))
        character(*), parameter :: this_routine = 'sort_arr_int'

        !/* The size of the array to sort */
        n = size(arr(1,:))

        nstack = 0
        lo = 1
        hi = n
        do while (.true.)
        !    /* If the section/partition we are looking at is smaller than
            !   nInsertSort then perform an insertion sort */
            if (hi - lo < nInsertionSort) then
                do j = lo + 1, hi
                    a = arr(:,j)
                    do i = j - 1, 1, -1
                        if (arr_lt(arr(:,i), a)) exit
                        arr(:,i+1) = arr(:,i)
                    enddo
                    arr(:,i+1) = a
                enddo

                if (nstack == 0) exit
                hi = stack(nstack)
                lo = stack(nstack-1)
                nstack = nstack - 2

        !    /* Otherwise start partitioning with quicksort. */
            else
                ! Pick a partitioning element, and arrange such that
                ! arr(lo) <= arr(lo+1) <= arr(hi)
                pivot = (lo + hi) / 2
                call swap (arr(:,pivot), arr(:,lo + 1))
                if (arr_gt(arr(:,lo), arr(:,hi))) &
                    call swap (arr(:,lo), arr(:,hi))
                if (arr_gt(arr(:,lo+1), arr(:,hi))) &
                    call swap (arr(:,lo+1), arr(:,hi))
                if (arr_gt(arr(:,lo), arr(:,lo+1))) &
                    call swap (arr(:,lo), arr(:,lo+1))

                i = lo + 1
                j = hi
                a = arr(:,lo + 1) ! a is the pivot value
                do while (.true.)
        !            /* Scand down list to find element > a */
                    i = i + 1
                    do while (arr_lt(arr(:,i), a))
                        i = i + 1
                    enddo

        !            /* Scan down list to find element < a */
                    j = j - 1
                    do while (arr_gt(arr(:,j), a))
                        j = j - 1
                    enddo

        !            /* When the pointers crossed, partitioning is complete. */
                    if (j < i) exit

        !            /* Swap the elements, so that all elements < a end up
        !               in lower indexed variables. */
                    call swap (arr(:,i), arr(:,j))
                enddo

                ! Insert partitioning element
                arr(:,lo + 1) = arr(:,j)
                arr(:,j) = a

        !        /* Push the larger of the partitioned sections onto the stack
        !           of sections to look at later.
        !           --> need fewest stack elements. */
                nstack = nstack + 2
                if (nstack > nStackMax) call stop_all (this_routine, &
                                        "parameter nStackMax too small")
                if (hi - i + 1 >= j - 1) then
                    stack(nstack) = hi
                    stack(nstack-1) = i
                    hi = j - 1
                else
                    stack(nstack) = j - 1
                    stack(nstack-1) = lo
                    lo = i
                endif
            endif
        enddo

        
    end subroutine

    subroutine sort_int (arr, n)

        ! Perform a quicksort on an array of integers, arr(n). Uses the 
        ! sample code in NumericalRecipies as a base. */

        integer, intent(in) :: n
        integer, intent(inout) :: arr(n)
        integer, parameter :: nInsertionSort = 7
        ! sort needs auxiliary storage of length 2*log_2(n). */
        integer, parameter :: nStackMax = 50
        character(*), parameter :: this_routine = 'sort_int'

        integer :: pivot, lo, hi, i, j, a
        ! Oh, how lovely it would be to be able to use push/pop and not need
        ! to guess a size of the stack to start with */
        integer :: stack(nStackMax), nstack


        nstack = 0
        lo = 1
        hi = n
        do while (.true.)
            ! If the section/partition we are looking at is smaller than
            ! nInsertSort then perform an insertion sort */
            if (hi - lo < nInsertionSort) then
                do j = lo + 1, hi
                    a = arr(j)
                    do i = j - 1, 1, -1
                        if (arr(i) < a) exit
                        arr(i+1) = arr(i)
                    enddo
                    arr(i+1) = a
                enddo

                if (nstack == 0) exit
                hi = stack(nstack)
                lo = stack(nstack-1)
                nstack = nstack - 2

            ! Otherwise start partitioning with quicksort. */
            else
                ! Pick a partitioning element, and arrange such that
                ! arr(lo) <= arr(lo+1) <= arr(hi)
                pivot = (lo + hi) / 2
                call swap (arr(pivot), arr(lo + 1))

                if (arr(lo) > arr(hi)) call swap (arr(lo), arr(hi))
                if (arr(lo+1) > arr(hi)) call swap (arr(lo+1), arr(hi))
                if (arr(lo) > arr(lo+1)) call swap (arr(lo), arr(lo+1))

                i = lo + 1
                j = hi
                a = arr(lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scand down list to find element > a */
                    i = i + 1
                    do while (arr(i) < a)
                        i = i + 1
                    enddo

                    ! Scan down list to find element < a */
                    j = j - 1
                    do while (arr(j) > a)
                        j = j - 1
                    enddo

                    ! When the pointers crossed, partitioning is complete. */
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables. */
                    call swap (arr(i), arr(j))
                enddo

                ! Insert partitioning element
                arr(lo + 1) = arr(j)
                arr(j) = a

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements. */
                nstack = nstack + 2
                if (nstack > nStackMax) call stop_all (this_routine, &
                                        "parameter nStackMax too small")
                if (hi - i + 1 >= j - 1) then
                    stack(nstack) = hi
                    stack(nstack-1) = i
                    hi = j - 1
                else
                    stack(nstack) = j - 1
                    stack(nstack-1) = lo
                    lo = i
                endif
            endif
        enddo

        
    end subroutine



end module

