program fruit_driver
    use fruit
    use example_test

    call init_fruit
    call example_test_driver
    call fruit_summary
    call fruit_finalize
end program fruit_driver
