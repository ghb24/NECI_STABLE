title: KP-FCIQMC Block 
---

### KP-FCIQMC Block

[TOC]

This block enables the Krylov-projected FCIQMC (KPFCIQMC) method [4]
which is fully implemented in NECI. It requires `dneci` or `mneci` to be
run. When specifying the KP-FCIQMC block, the METHODS block should be
omitted. This block is started with the `kp-fciqmc` keyword and
terminated with the `end-kp-fciqmc` keyword.

-   **kp-fciqmc**  
    Starts the KP-FCIQMC block

-   **end-kp-fciqmc**  
    Terminates the KP-FCIQMC block.

-   **num-krylov-vecs \(N\)**  
    \(N\) specifies the total number of Krylov vectors to sample.

-   **num-iters-between-vecs \(N\)**  
    \(N\) specifies the (constant) number of iterations between each
    Krylov vector sampled. The first Krylov vector is always the
    starting wave function.

-   **num-iters-between-vecs-vary \(i_{12}\), \(i_{23}\),
    \(i_{34}\)...**  
    \(i_{n,n+1}\) specifies the number of iterations between the nth and
    (n+1)th Krylov vectors. The number of parameters input should be the
    number of Krylov vectors asked for minus one. The first Krylov
    vector is always the starting wave function.

-   **num-repeats-per-init-config \(N\)**  
    \(N\) specifies the number repeats to perform for each initial
    configuration, i.e. the number of repeats of the whole evolution,
    from the first sampled Krylov vector to the last. The projected
    Hamiltonian and overlap matrix estimates will be output for each
    repeat, and the averaged values of these matrices used to compute
    the final results.

-   **averagemcexcits-hamil \(N\)**  
    When calculating the projected Hamiltonian estimate, an FCIQMC-like
    spawning is used, rather than calculating the elements exactly,
    which would be too computationally expensive. Here, \(N\) specifies
    the number of spawnings to perform from each walker from each Krylov
    vector when calculating this estimate. Thus, increasing \(N\) should
    improve the quality of the Hamiltonian estimate.

-   **finite-temperature**  
    If this option is included then a finite-temperature calculation is
    performed. This involves starting from several different random
    configurations, whereby walkers are distributed on random
    determinants. The number of initial configurations should be
    specified with the num-init-configs option.

-   **num-init-configs \(N\)**  
    \(N\) specifies the number of initial configurations to perform the
    sampling over. An entire FCIQMC calculation will be performed, and
    an entire subspace generated, for each of these configurations. This
    option should be used with the finite-temperature option, but is not
    necessary for spectral calculations where one always starts from the
    same initial vector.

-   **memory-factor \(x\)**  
    This option is used to specify the size of the array allocated for
    storing the Krylov vectors. The number of slots allocated to store
    unique determinants in the array holding all Krylov vectors will be
    equal to \(ABx\), where here \(A\) is the length of the main walker
    list, \(B\) is the number of Krylov vectors, and \(x\) is the value
    input with this option.

-   **num-walker-per-site-init \(x\)**  
    For finite-temperature jobs, \(x\) specifies the number of walkers
    to place on a determinant when it is chosen to be occupied.

-   **exact-hamil**  
    If this option is specified then the projected Hamiltonian will be
    calculated exactly for each set of Krylov vectors sampled, rather
    than randomly sampling the elements via an FCIQMC-like spawning
    dynamic.

-   **fully-stochastic-hamil**  
    If this option is specified then the projected Hamiltonian will be
    estimated without using the semi-stochastic adaptation. This will
    decrease the quality of the estimate, but may be useful for
    debugging or analysis of the method.

-   **init-correct-walker-pop**  
    For finite-temperature calculations on multiple cores, the initial
    population may not be quite as requested. This is because the
    quickest (and default) method involves generating determinants
    randomly and sending them to the correct processor at the end. It is
    possible in this process that walkers will die in annihilation.
    However, if this option is specified then each processor will throw
    away spawns to other processors, thus allowing the correct total
    number of walkers to be spawned.

-   **init-config-seeds seed1, seed2...**  
    If this option is used then, for finite-temperature calculations, at
    the start of each calculation over an initial configuration, the
    random number generator will be re-initialised with the
    corresponding input seed. The number of seeds provided should be
    equal to the number of initial configurations.

-   **all-sym-sectors**  
    If this option is specified then the FCIQMC calculation will be run
    in all symmetry sectors simultaneously. This is an option relevant
    for finite-temperature calculations.

-   **scale-population**  
    If this option is specified then the initial population will be
    scaled to the populaton specified with the ‘totalwalkers’ option in
    the Calc block. This is relevant for spectral calculations when
    starting from a perturbed `POPSFILE` wave function, where the
    initial population is not easily controlled.

In spectral calculations, one also typically wants to consider a
particular perturbation operator acting on the ground state wave
functions. Therefore, you must first perform an FCIQMC calculation to
evolve to the ground state and output a `POPSFILE`. You should then
start the KP-FCIQMC calculation from that `POPSFILE`. To apply a
perturbation operator to the `POPSFILE` wave function as it is read in,
use the pops-creation and pops-annihilate options. These allow operators
such as
\[\hat{V} = \hat{c}_i \hat{c}_j + \hat{c}_k \hat{c}_l \hat{c}_m\] to be
applied to the `POPSFILE` wave function. The general form is
`pops-annihilate n_sum` `orb1 orb2...` `...` where \(n_{sum}\) is the
number of terms in the sum for \(\hat{V}\) (2 in the above example), and
\(orbi\) specify the spin orbital labels to apply. The number of lines
of such orbitals provided should be equal to \(n_{sum}\). The first line
provides the orbital labels for the first term in the sum, the second
line for the second term, etc...

