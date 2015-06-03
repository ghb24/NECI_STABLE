program fruit_driver
    use fruit
    use example_test

    use det_bit_ops_tests

    call init_fruit
    call example_test_driver
    call det_bit_ops_drive_tests
    call fruit_summary
    call fruit_finalize
end program fruit_driver
