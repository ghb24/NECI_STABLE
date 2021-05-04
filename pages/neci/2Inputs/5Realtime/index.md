title: Realtime Block 
---
### REALTIME Block

The REALTIME block enables the calculation of the real-time evolution of
a given state and is only required for real-time calculations. Real-time
evolution strictly requires `kneci` or `kmneci` and using `kmneci` is
strongly recommended. This block is started with the `realtime` keyword
and terminatd with the `endrealtime` keyword.

-   **realtime**  
    Starts the realtime block. This automatically enables `readpops`,
    and reads from a file named `<popsfilename>.0`.

-   **endrealtime**  
    Terminates the REALTIME block.

-   **single \(i\) \(a\)**  
    Applies a single excitation operator exciting from orbital \(i\) to
    orbital \(a\) to the initial state before starting the calculation.

-   **lesser \(i\) \(j\)**  
    Calculates the one-particle Green’s function \(<c_i^\dagger(t)
          c_j>\). Requires using one electron less than in the popsfile.

-   **greater \(i\) \(j\)**  
    Calculates the one-particle Green’s function \(<c_i(t)
          c_j^\dagger>\). Requires using one electron more than in the
    popsfile.

-   **start-hf**  
    Do not read in a popsfile but start from a single determinant.

-   **rotate-time \(\alpha\)**  
    Calculates the evolution along a trajectory \(e^{i\alpha}t\) instead
    of pure real-time trajectory.

-   **dynamic-rotation [\(\eta\)]**  
    Determine a time-dependent \(\alpha\) on the fly for `rotate-time`.
    \(\eta\) is optional and a damping parameter, defaults to 0.05.
    \(\alpha\) is then chosen such that the walker number remains
    constant.

-   **rotation-threshold \(N\)**  
    Grow the population to \(N\) walkers before starting to adjust
    \(\alpha\).

-   **stepsalpha \(n\)**  
    The number of timesteps between two updates of the \(\alpha\)
    parameter when using `dynamic-rotation`.

-   **log-trajectory**  
    Output the time trajectory to a separate file. This is useful if the
    same calculation shall be reproduced.

-   **read-trajectory**  
    Read the time trajectory from disk, using a file created by NECI
    with the `log-trajectory` keyword.

-   **live-trajectory**  
    Read the time trajectory from disk, using a file which is currently
    being created by NECI with the `log-trajectory` keyword. Can be used
    to create additional data for the same trajectory while the original
    calculation is still running.

-   **noshift**  
    Do not apply a shift during the real-time evolution. Strongly
    recommended.

-   **stabilize-walkers [\(S\)]**  
    Use the shift to stabilize the walker number if it drops below 80%
    of the peak value. \(S\) Is optional and is an asymptotic value used
    to fix the shift.

-   **energy-benchmark \(E\)**  
    Set the energy origin to \(E\) by applying a global, constant shift
    to the Hamiltonian. Can be chosen arbitrarily, but a reasonable
    selection can greatly help efficiency.

-   **rt-pops**  
    A second popsfile is supplied containing a time evolved state
    created with the `realtime` keyword whose time evolution is to be
    continued. In this case, the original popsfile is still required for
    calculating properties, so two popsfiles will be read in.
