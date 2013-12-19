module stoch_lanczos

    use stoch_lanczos_procs

    implicit none

contains

    subroutine perform_stoch_lanczos(lanczos)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer :: iconfig, irun, ivec, iiter

        call init_stoch_lanczos(lanczos)

        do iconfig = 1, lanczos%nconfigs

            do irun = 1, lanczos%nruns

                call create_initial_config(lanczos, irun)

                do ivec = 1, lanczos%nkrylov_vecs
                    do iiter = 1, lanczos%niters
                    end do
                end do
            end do
        end do

    end subroutine perform_stoch_lanczos

end module stoch_lanczos
