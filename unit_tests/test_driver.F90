program fruit_driver
    use fruit

    use det_bit_ops_tests

    integer :: failed_count

    call init_fruit
    call det_bit_ops_drive_tests
    call fruit_summary
    call fruit_finalize

    ! If we have failed any tests, we ought to return return a failing return code
    ! --> This should get picked up by Jenkins
    call get_failed_count(failed_count)
    if (failed_count /= 0) then
        stop -1
    end if
end program fruit_driver
