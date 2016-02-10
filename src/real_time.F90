#include "macros.h"

module real_time

    use real_time_init, only: init_real_time_calc_single

    implicit none

! main module file for the real-time implementation of the FCIQMC algorithm
! created on 04.02.2016 by Werner Dobrautz

! first i have to do a brainstorm on how to implement Ole's idea of the 
! real time version and also discuss with Ali, Simon and George to avoid 
! unnecessary effort

! implementation idea: 
! sample the real time Schroedinger equation by integrating:
! i d/dt |Psi(t)> = (H - E0 - ie)|Psi(t)> 

! with the 2nd order runge-kutta method: 
! d/dt y(t) = f(t,y)

! -> y(n+1) = y(n) + h f(n dt, y(n))
! t = n dt
! y(n+1) = y(n) + k2
! k1 = dt f(n dt, y(n)) 
! k2 = dt f(n + dt/2, y(n) + k1/2) 

! what does that mean in the dynamics? 

! from brainstorming on the 2nd order runge kutta implementation of the 
! original imaginary time FCIQMC algorithm, a few clarifications surfaced: 

! the RK2 needs the application of the square of the Hamiltonian essentially
! this implies huge changes to the underlying machinery in the current 
! FCIQMC implementation, which often relies on the type of exciations
! (single, doubles) 
! since the intermediatly created list y(n) + k1 / 2 already contains 
! single and double excitations, spawning from this list, increases the 
! possible excitations in the final y(n) + k2 -> y(n+1) to 
! singles, double, triples and quadrupels

! but the essential goal in the real-time code will be to obtain the 
! overlap: 
! <y(0)|y(t)> 

! to obtain the spectral information of the system up to a certain time 
! t_max

! i should devise an action plan to optimally implement this method 

! first of all i definetly need to use kneci to be able to handle the 
! imaginary and real parts of the equation
! dy(t)/dt = -i(H - E0 -ie)y(t) 

! then, since Ole showed already, and what i already know, back from the 
! work in Graz, a averaging over multiple runs, starting from different 
! stochastic representations of the converged ground state wave function is 
! necessary! 

! so probably a multiple kneci version, as developed by Dongxia is needed

! so whats the ideas?

! we have to start from a converged imaginary-time ground state wave-function! 
! a way could be to use a printed result from a previous calculationm through
! a popsfile. 
! but since we want average over multiple calculations, this would need 
! unnecessary amount of storage, I/O etc.
! but we should definetly be able to have this option! in the case the 
! imaginary-time calculation is very involved and already takes a lot of time 
! so maybe be able to read in #n popsfiles, printed at different times of 
! the imag-time calculation and start that amount of mkneci processes, which 
! do the real-time calculation then. 

! or run a "normal" imag-time run and then after convergence, start the 
! #n real-time calculations from diffferent starting points, as indicated by 
! Ole in this unterlagen
! or start with running multiple #n mneci runs with different seeds for the 
! random number generator, and then at the same time in equilibration, start 
! the real-time calculation

! this would imply a completly seperated code-basis in which the real-time 
! implementation is performed. which probably is a good idea anyway to keep 
! the code clean. since we need totally different information and stats 
! anyway..

! at each time step, we could then average <y(0)|y(t)> from the different 
! mkneci processes and just store one quantity and additional statistical info

! back in Graz i also came to the conclusion, that calculating the overlap 
! between different runs was of great help -> we could do that here too

! on the implementation of the new dynamics: 
! the underlying differential equation changes to (for the greater GF)
! dy(t)/dt = -i(H - E0 - ie)y(t) 
! this means 
! 1) we need both real and imaginary walkers 
! 2) when using the 2nd order RK method: 
!       y(n+1) = y(n) + k2(n)
!       k1(n) = -i dt (H - E0 - ie)y(n)
!       k2(n) = -i dt (H - E0 - ie)(y(n) + k1 / 2) 
!       
!       we could implement that by creating an intermediate determinant 
!       list obtained by the spawn(+annihilation) and death/cloning step
!       from y(n) + k1 / 2 
!       and based upon the spawns from that list we can combine that to the 
!       new y(n+1) = y(n) + k2 for the list at the next time-step
! 
!   or we can work out the final recursion formula:
!       y(n+1) = [(1 - e dt)(1 - i dt H') - dt^2 / 2 H'^2] y(n) (see doc.)
!
!   and do the spawning/cloning/death etc. directly based on this equation
!   lets see and talk to ali/ole about that 

! so the workflow of real-time propagations will be: 

! 1) run a normal imaginary time propagation until convergence and print out 
!       #n popsfiles as the groundstate information

! 2) terminate the imaginary run and start the real-time propagation with #n
!       threads mneci style

! 3) specify a set of specific orbitals j, on which either a^+_j or a_j acts
!       with possible symmetry constraints also, or specify k-space calculation
!       like for the Hubbard or UEG model. 
!       for each j or k also the coresponding <0|a_i/a^+_i or multiple if 
!       possible act, to calculate the overlaps. applying multiple <0|a_i is 
!       cheap, so all of them should be done. (for k-space k'=k due to symmetry)

! 4) do the actual real time propagation! 

contains

    subroutine perform_real_time_fciqmc
        ! main real-time calculation routine
        ! do all the setup, read-in and calling of the "new" real-time MC loop
        implicit none

        character(*), parameter :: this_routine = "perform_real_time_fciqmc"

        ! call the real-time setup routine and all the initialization
        call init_real_time_calc_single()

        ! rewrite the major original neci core loop here and adapt it to 
        ! the new necessary real-time stuff
        ! check nicks kp code, to have a guideline in how to go into that! 

    end subroutine perform_real_time_fciqmc
end module real_time

! wrapper (dont know why this is necessary quite..)
subroutine perform_real_time_fciqcm_wrap
    use real_time, only: perform_real_time_fciqmc
    implicit none

    call perform_real_time_fciqmc

end subroutine perform_real_time_fciqcm_wrap
