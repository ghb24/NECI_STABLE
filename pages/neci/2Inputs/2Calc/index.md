title: CALC Block
---

### CALC Block

[TOC]

The CALC block is used to set options concerning the simulation
parameters and modes of FCIQMC. The block starts with the `calc` keyword
and ends with the `endcalc` keyword.

-   <span style="color: red">**calc**  
    </span> Starts the CALC block

-   <span style="color: red">**endcalc**  
    </span> Terminates the CALC block

-   <span style="color: blue">**time \(t\)**  
    </span> Set the maximum time \(t\) in minutes the calculation is
    allowed to run. After \(t\) minutes, the calculation will end.

-   <span style="color: blue">**nmcyc \(n\)**  
    </span> Set the maximum number of iterations the calculation is
    allowed to do. After \(n\) iterations, the calculation will end.

-   **eq-cyc \(n\)**  
    Set the number of iterations to perform after reaching variable
    shift mode. \(n\) Iterations after entering the variable shift mode,
    the calculation will end. If both `eq-cyc` and `nmcyc` are given,
    the calculation ends when either of the iteration limits is reached.

-   **seed \(s\)**  
    Sets the seed of the random number generator to \(s\). This can be
    used to specifically probe for stochastic effects, but is generally
    not required.

-   **averageMcExcits \(x\)**  
    Sets the average number of spawning attempts from each walker to
    \(x\).

-   **rdmSamplingIters \(n\)**  
    Set the maximum number of iterations used for sampling the RDMs to
    \(n\). After \(n\) iterations of sampling RDMs, the calculation will
    end.

-   **load-balance-blocks [OFF]**  
    Distribute the determinants blockwise in a dynamic fashion to
    maintain equal load for all processors. This is enabled by default
    and has one optional argument OFF. If given, the load-balancing is
    disabled.

-   **energy**  
    Additionally calculate and print the ground state energy using an
    exact diagonalization technique.

-   **averageMcExcits \(n\)**  
    The number of spawns to attempt per walker. Defaults to \(1\) and
    should not be changed without good reason.

-   **adjust-averageMcExcits**  
    Dynamically update the number of spawns attempted per walker. Can be
    used if the excitation generator creates a lot of invalid
    excitations, but should be avoided else.

-   **scale-spawns**  
    Store the maximum value of \(\frac{H_{ij}}{p_{gen}}\) for each
    determinant and use it to estimate the number of spawns per walker
    to prevent blooms. Useful when this fraction strongly depends on the
    determinant.

#### Population control options

-   <span style="color: red">**totalWalkers \(n\)**  
    </span> Sets the targeted number of walkers to \(n\). This means,
    the shift will be varied to keep the walker number constant once it
    reaches \(n\).

-   <span style="color: blue">**diagShift \(S\)**  
    </span> Set the initial value of the shift to \(S\). A value of
    \(S<0\) is not recommended, as it will decrease the population from
    the beginning.

-   <span style="color: blue">**shiftDamp \(\zeta\)**  
    </span> Set the damping factor used in the shift update scheme to
    \(\zeta\). Defaults to \(10\).

-   <span style="color: blue">**stepsSft \(n\)**  
    </span> Sets the number of steps per update cycle of the shift to
    \(n\). Defaults to \(100\).

-   **fixed-n0 \(n_0\)**  
    Instead of varying the shift to fix the total number of walkers,
    keep the number of walkers at the reference fixed at \(n_0\).
    Automatically sets `stepsSft 1` and overwrites any `stepssft`
    options given.

-   **targetGrowRate \(grow\) \(walks\)**  
    When the number of walkers in the calculation exceeds \(walk\), the
    shift is iteratively adjusted to maintain a fixed grow rate \(grow\)
    until reaching the requested number of total walkers.

-   **jump-shift [OFF]**  
    When entering the variable shift mode, the shift will be set to the
    current projected energy. This is enabled by default. There is an
    optional argument OFF that disables this behaviour.

-   **pops-jump-shift**  
    Reset the shift when restarting a previous calculation to the
    current projected energy instead of using the shift from the
    previous calculation.

-   **trunc-nopen \(n\)**  
    Restrict the Hilbert space of the calculation to those determinants
    with at most \(n\) unpaired electrons.

-   **avGrowthRate [OFF]**  
    Average the change in walker number used to calculate the shift.
    This is enabled by default and has one optional argument OFF, which,
    when given, turns the option off.

#### Real walker coefficient options

-   <span style="color: blue">**allRealCoeff**  
    </span> Allow determinants to have non-integer population. There is
    a minimal population below which the population of a determinant
    will be rounded stochastically. This defaults to \(1\).

-   <span style="color: blue">**realSpawnCutoff \(x\)**  
    </span> Continuous real spawning will be performed, unless the spawn
    has weight less than x. In this case, the weight of the spawning
    will be stochastically rounded up to x or down to zero, such that
    the average weight of the spawning does not change. This is a method
    of removing very low weighted spawnings from the spawned list, which
    require extra memory, processing and communication. A reasonable
    value for x is 0.01.

-   **realCoeffbyExcitLevel \(n\)**  
    Allow all determinants up to an excitation level of \(n\) to have
    non-integer population.

-   **setOccupiedThresh \(x\)**  
    Set the value for the minimum walker weight in the main walker list.
    If, after all annihilation has been performed, any determinants have
    a total weight of less than x, then the weight will be
    stochastically rounded up to x or down to zero such that the average
    weight is unchanged. Defaults to \(1\), which should only be changed
    with good reason.

-   **energy-scaled-walkers [\(mode\) \(\alpha\) \(\beta\)]**  
    Scales the occupied threshold with the energy of a determinant. Has
    three optional arguments and requires `allRealCoeff`. The argument
    \(mode\) can be one of EXPONENTIAL, POWER, EXP-BOUND or NEGATIVE and
    defaults to POWER. Both \(\alpha\) and \(\beta\) default to 1 and
    `realspawncutoff` \(\beta\) is implied.

#### Time-step options

-   **<span style="color: red">tau \(\tau\)</span> [SEARCH]**  
    Sets the timestep per iteration to \(\tau\). Has one optional
    argument SEARCH. If given, the time-step will be iteratively updated
    to keep the calculation stable.

-   <span style="color: blue">**hist-tau-search \[\(c\) \(nbins\)
    \(bound\)\]**  
    </span> Update the time-step based on histogramming of the ratio
    \(\frac{H_{ij}}{p(i|j)}\). Not compatible with the tau \(\tau\)
    SEARCH option. The three arguments \(c\), \(nbins\) and \(bound\)
    are optional. \(0<c<1\) is the fraction of the histogram used for
    determining the new timestep, \(nbins\) the number of bins in the
    histogram and \(bound\) is the maximum value of
    \(\frac{H_{ij}}{p(i|j)}\) to be stored.  
    For spin-adapted GUGA calculations this option is *highly*
    recommended! Otherwise the time-step can become quite small in these
    simulations.

-   <span style="color: blue">**max-tau \(\tau_\text{max}\)**  
    </span> Sets the maximal value of the time-step to
    \(\tau_\text{max}\). Defaults to \(1\).

-   **min-tau [\(\tau_\text{min}\)]**  
    Sets the minimal value of the time-step to \(\tau_\text{min}\) and
    enables the iterative update of the time-step. Defaults to
    \(10^{-7}\). The argument \(\tau_\text{min}\) is optional.

-   **keepTauFixed**  
    Do never update \(\tau\) and the related parameter
    \(p_\text{singles}\), \(p_\text{doubles}\) or \(p_\text{parallel}\).

-   **truncate-spawns [\(n\) UNOCC]**  
    Truncate spawns which are larger than a threshold value \(n\). Both
    arguments are optional, \(n\) defaults to \(3\). If UNOCC is given
    the truncation is restricted to spawns onto unoccupied. Useful in
    combination with hist-tau-search.

-   **maxWalkerBloom \(n\)**  
    The time step is scaled such that at most \(n\) walkers are spawned
    in a single attempt, with the scaling being guessed from previous
    spawning attempts.

#### Wave function initialization options

-   <span style="color: blue">**walkContGrow**  
    </span> When reading in a wave function from a file, do not set the
    shift or enter variable shift mode.

-   <span style="color: blue">**defineDet \(det\)**  
    </span> Sets the reference determinant of the calculation to
    \(det\). If no other initialisation is specified, this will also be
    the initial wave function. The format can either be a
    comma-separated list of spin orbitals, a range of spin orbitals
    (like 12-24) or a combination of both.  
    For spin-adapted calculations using the `GUGA` keyword defining a
    starting reference CSF manually is *highly* encouraged, as the
    automatic way often fails. It works in similar ways as for SDs,
    however odd numbered singly occupied orbitals indicate a positive
    spin coupled orbital in GUGA CSFs and even numbered singly occupied
    orbitals negatively spin coupled. So one needs to be careful to not
    define a unphysical CSF with an negative intermediate total spin
    \(S_i < 0\). E.g. a CSF like:  
    `definedet 1 2 4 5`  
    would cause the calculation to crash, as the first singly occupied
    orbital (4) would cause the total spin to be negative \(S_i < 0\).
    When defining the starting CSF one also needs to ensure that the
    defined CSFs satisfies the total number of electrons and total \(S\)
    defined in the `System` block of the input with the keywords `nel`
    and `guga`. As an example a triplet (`guga 2`) CSF with 4 electrons
    (`nel 4`) would be  
    `definedet 1 2 3 5`  
    with the first spatial orbital doubly occupied and 2 open-shell
    orbitals (positively spin coupled, hence odd numbers).

-   **readPops**  
    Read in the wave function from a file and use the read-in wave
    function for initialisation. In addition to the wave function, also
    the time-step and the shift are read in from the file. This starts
    the calculation in variable shift mode, maintaining a constant
    walker number, unless `walkContGrow` is given.

-   **readpops-changeref**  
    Allow the reference determinant to be updated after reading in a
    wave function.

-   **startSinglePart [\(n\)]**  
    Initialise the wave function with \(n\) walkers on the reference
    only unless specified differently. The argument \(n\) is optional
    and defaults to \(1\).

-   **proje-changeRef [\(frac\) \(min\)]**  
    Allow the reference to change if a determinant obtains \(frac\)
    times the population of the current reference and the latter has a
    population of at least \(min\). Both arguments are optional and
    default to \(1.2\) and \(50\), respectively. This is enabled by
    default.

-   **no-changeref**  
    Never change the reference.

#### Initiator options

-   <span style="color: blue">**truncInitiator**  
    </span> Use the initiator method [3].

-   <span style="color: blue">**addToInitiator \(x\)**  
    </span> Sets the initiator threshold to \(x\), so any determinant
    with more than \(x\) walkers will be an initiator.

-   **senior-initiators [\(age\)]**  
    Makes any determinant that has a half-time of at least \(age\)
    iterations an initiator. \(age\) is optional and defaults to \(1\).

-   **superInitiator [\(n\)]**  
    Create a list of \(n\) superinitiators, from which all connected
    determinants are set to be initiators. The superinitiators are
    chosen according to population. \(n\) is optional and defaults to
    \(1\).

-   **coherent-superInitiators [\(mode\)]**  
    Apply a restriction on the sign-coherence between a determinant and
    any connected superinitiator to determine whether it becomes an
    initiator due to connection. The optional argument \(mode\) can be
    chosen from STRICT, WEAK, XI, AV and OFF. The default is WEAK and is
    enabled by default if the `superInitiators` keyword is given.

-   **dynamic-superInitiators \(n\)**  
    Updates the list of superinitiators every \(n\) steps. This is
    enabled by default with \(n=100\) if the `superInitiators` keyword
    is given. A value of 0 indicates no update. This implies the
    `dynamic-core n` option (with the same \(n\)) unless specified
    otherwise.

-   **allow-signed-spawns \(mode\)**  
    Never abort spawns with a given sign, regardless of initiators.
    \(mode\) can be either POS or NEG, indicating the sign to keep.

-   **initiator-space**  
    Define all determinants within a given initiator space as
    initiators. The space is specified through one of the following
    keywords

    -   **doubles-initiator**  
        Use the reference determinant and all single and double
        excitations from it to form the initiator space.

    -   **cas-initiator cas1 cas2**  
        Use a CAS space to form the initiator space. The parameter cas1
        specifies the number of electrons in the cas space and cas2
        specifies the number of virtual spin orbitals (the cas2 highest
        energy orbitals will be virtuals).

    -   **ras-initiator ras1 ras2 ras3 ras4 ras5**  
        Use a RAS space to form the initiator space. Suppoose the list
        of spatial orbitals are split into three sets, RAS1, RAS2 and
        RAS 3, ordered by their energy. ras1, ras2 and ras3 then specify
        the number of spatial orbitals in RAS1, RAS2 and RAS3. ras4
        specifies the minimum number of electrons in RAS1 orbitals. ras5
        specifies the maximum number of electrons in RAS3 orbitals.
        These together define the RAS space used to form the initiator
        space.

    -   **optimised-initiator**  
        Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
        230201). One also needs to use either
        optimised-initiator-cutoff-amp or optimised-initiator-cutoff-num
        with this option.

    -   **optimised-initiator-cutoff-amp \(x1\), \(x2\), \(x3\)...**  
        Perform the optimised initiator option, and in iteration \(i\),
        choose which determinants to keep by choosing all determinants
        with an amplitude greater than \(xi\) in the ground state of the
        space (see PRL 109, 230201). The number of iterations is
        determined by the number of parameters provided.

    -   **optimised-initiator-cutoff-num \(n1\), \(n2\), \(n3\)...**  
        Perform the optimised initiator option, and in iteration \(i\),
        choose which determinants to keep by choosing the ni most
        significant determinants in the ground state of the space (see
        PRL 109, 230201). The number of iterations is determined by the
        number of parameters provided.

    -   **fci-initiator**  
        Use all determinants to form the initiator space. A fully
        deterministic projection is therefore performed with this
        option.

    -   **pops-initiator \(n\)**  
        When starting from a POPSFILE, this option will use the \(n\)
        most populated determinants from the popsfile to form the
        initiator space.

    -   **read-initiator**  
        Use the determinants in the INITIATORSPACE file to form the
        initiator space. A INITIATORSPACE file can be created by using
        the write-initiator option in the LOGGING block.

#### Adaptive shift options

-   **auto-adpative-shift [\(t\) \(\alpha\) \(c\)]**  
    Scale the shift per determinant based on the acceptance rate on a
    determinant. Has three optional arguments. The first is the
    threshold value \(t\) which is the minimal number of spawning
    attempts from a determinant over the full calculation required
    before the shift is scaled, with a default of \(10\). The second is
    the scaling exponent \(\alpha\) with a default of \(1\) and the last
    is the minimal scaling factor, which uses
    \(\frac{1}{\text{HF conn.}}\) as default.

-   **linear-adaptive-shift [\(\sigma\) \(f_1\) \(f_2\)]**  
    Scale the shift per determinant linearly with the population of a
    determinant. All arguments are optional and define the function used
    for scaling. \(\sigma\) gives the minimal walker number required to
    have a shift and defaults to \(1\), \(f_1\) the shift fraction to be
    applied at \(\sigma\) with a default of \(0\) and \(f_2\) is the
    shift fraction to be applied at the initiator threshold, defaults to
    \(1\). Every initiator is applied the full shift.

-   **exp-adaptive-shift [\(\alpha\)]**  
    Scales the shift expoentially with the population of a determinant.
    The optional argument \(\alpha\) is the exponent of scaling, the
    default is \(2\).

-   **core-adaptive-shift**  
    By default, determinants in the corespace are always applied the
    full shift. Using this option also scales the shift in the
    corespace.

-   **aas-matele2**  
    Uses the matrix elements for determining the scaling factor in the
    `auto-adaptive-shift`. The recommended option to scale the shift.

#### Multi-replica options

-   **multiple-initial-refs**  
    Define a reference determinant for each replica. The following \(n\)
    lines give the reference determinants as comma-separated lists of
    orbitals, where \(n\) is the number of replicas.

-   **orthogonalise-replicas**  
    Orthogonalise the replicas after each iteration using Gram Schmidt
    orthogonalisation. This will converge each replica to another state
    in a set of orthogonal eigenstates. Can be used for excited state
    search.

-   **orthogonalise-replicas-symmetric**  
    Use the symmetric Löwdin orthonaliser instead of Gram Schmidt for
    orthogonalising the replicas.

-   **replica-single-det-start**  
    Starts each replica from a different excited determinant.

#### Semi-stochastic options

-   <span style="color: blue">**semi-stochastic**  
    </span> Turn on the semi-stochastic adaptation.

-   <span style="color: blue">**pops-core \(n\)**  
    </span> This option will use the \(n\) most populated determinants
    to form the core space.

-   **doubles-core**  
    Use the reference determinant and all single and double excitations
    from it to form the core space.

-   **cas-core cas1 cas2**  
    Use a CAS space to form the core space. The parameter cas1 specifies
    the number of electrons in the cas space and cas2 specifies the
    number of virtual spin orbitals (the cas2 highest energy orbitals
    will be virtuals).

-   **ras-core ras1 ras2 ras3 ras4 ras5**  
    Use a RAS space to form the core space. Suppoose the list of spatial
    orbitals are split into three sets, RAS1, RAS2 and RAS 3, ordered by
    their energy. ras1, ras2 and ras3 then specify the number of spatial
    orbitals in RAS1, RAS2 and RAS3. ras4 specifies the minimum number
    of electrons in RAS1 orbitals. ras5 specifies the maximum number of
    electrons in RAS3 orbitals. These together define the RAS space used
    to form the core space.

-   **optimised-core**  
    Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
    230201). One also needs to use either optimised-core-cutoff-amp or
    optimised-core-cutoff-num with this option.

-   **optimised-core-cutoff-amp \(x1\), \(x2\), \(x3\)...**  
    Perform the optimised core option, and in iteration \(i\), choose
    which determinants to keep by choosing all determinants with an
    amplitude greater than \(xi\) in the ground state of the space (see
    PRL 109, 230201). The number of iterations is determined by the
    number of parameters provided.

-   **optimised-core-cutoff-num \(n1\), \(n2\), \(n3\)...**  
    Perform the optimised core option, and in iteration \(i\), choose
    which determinants to keep by choosing the ni most significant
    determinants in the ground state of the space (see PRL 109, 230201).
    The number of iterations is determined by the number of parameters
    provided.

-   **fci-core**  
    Use all determinants to form the core space. A fully deterministic
    projection is therefore performed with this option.

-   **read-core**  
    Use the determinants in the CORESPACE file to form the core space. A
    CORESPACE file can be created by using the write-core option in the
    LOGGING block.

-   **dynamic-core \(n\)**  
    Update the core space every \(n\) iterations, where \(n\) is
    optional and defaults to \(400\). This is enabled by default if the
    `superinitiators` option is given.

#### Trial wave function options

-   <span style="color: blue">**trial-wavefunction [\(n\)]**  
    </span> Use a trial wave function to obtain an estimate for the
    energy, as described in
    <a href="#sec:trial" data-reference-type="ref" data-reference="sec:trial">1.4</a>.
    The argument \(n\) is optional, when given, the trial wave function
    will be initialised \(n\) iterations after the variable shift mode
    started, else, at the start of the calculation. The trial wave
    function is defined through one of the following keywords

    -   <span style="color: blue">**pops-trial \(n\)**  
        </span> When starting from a POPSFILE, this option will use the
        \(n\) most populated determinants from the popsfile to form the
        trial space.

    -   **doubles-trial**  
        Use the reference determinant and all single and double
        excitations from it to form the trial space.

    -   **cas-trial cas1 cas2**  
        Use a CAS space to form the trial space. The parameter cas1
        specifies the number of electrons in the cas space and cas2
        specifies the number of virtual spin orbitals (the cas2 highest
        energy orbitals will be virtuals).

    -   **ras-trial ras1 ras2 ras3 ras4 ras5**  
        Use a RAS space to form the trial space. Suppoose the list of
        spatial orbitals are split into three sets, RAS1, RAS2 and RAS
        3, ordered by their energy. ras1, ras2 and ras3 then specify the
        number of spatial orbitals in RAS1, RAS2 and RAS3. ras4
        specifies the minimum number of electrons in RAS1 orbitals. ras5
        specifies the maximum number of electrons in RAS3 orbitals.
        These together define the RAS space used to form the trial
        space.

    -   **optimised-trial**  
        Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
        230201). One also needs to use either optimised-trial-cutoff-amp
        or optimised-trial-cutoff-num with this option.

    -   **optimised-trial-cutoff-amp \(x1\), \(x2\), \(x3\)...**  
        Perform the optimised trial option, and in iteration \(i\),
        choose which determinants to keep by choosing all determinants
        with an amplitude greater than \(xi\) in the ground state of the
        space (see PRL 109, 230201). The number of iterations is
        determined by the number of parameters provided.

    -   **optimised-trial-cutoff-num \(n1\), \(n2\), \(n3\)...**  
        Perform the optimised trial option, and in iteration \(i\),
        choose which determinants to keep by choosing the ni most
        significant determinants in the ground state of the space (see
        PRL 109, 230201). The number of iterations is determined by the
        number of parameters provided.

    -   **fci-trial**  
        Use all determinants to form the trial space. A fully
        deterministic projection is therefore performed with this
        option.

#### Memory options

-   **memoryFacPart \(x\)**  
    Sets the factor between the allocated space for the wave function
    and the required memory for the specified number of walkers to
    \(x\). Defaults to \(10\).

-   **memoryFacSpawn \(x\)**  
    Sets the factor between the allocated space for new spawns and the
    estimate of required memory for the spawns of the specified number
    of walkers on a single processor to \(x\). The memory required for
    spawns increases, the more processors are used, so when running with
    few walkers on relatively many processors, a large factor might be
    needed. Defaults to \(3\).

-   **prone-walkers**  
    Instead of terminating when running out of memory, randomly delete
    determinants with low population and few spawns.

-   **store-dets**  
    Employ extra memory to store additional information on the
    determinants that had to be computed on the fly else. Trades in
    memory for faster iterations.

#### Reduced density matrix (RDM) options

-   **rdmSamplingIters \(n\)**  
    Set the number of iterations for sampling the RDMs to \(n\). After
    \(n\) iterations of sampling, the calculation ends.

-   **inits-rdm**  
    Only take into account initiators when calculating RDMs. By default,
    only restricts to initiators for the right vector used in RDM
    calculation. This makes the RDMs non-variational, and the resulting
    energy is the projected energy on the initiator space.

-   **strict-inits-rdm**  
    Require both sides of the inits-rdm to be initiators.

-   **no-lagrangian-rdms**  
    This option disables the correction used for RDM calculation for the
    adaptive shift. Use this only for debugging purposes, as the
    resulting RDMs are flawed.

#### METHODS Block

The METHODS block is a subblock of CALC, i.e. it is specified inside the
CALC block. It sets the main algorithm to be used in the calculation.
The subblock is started with the `methods` keyword and terminated with
the `endmethods` keyword.

-   <span style="color: red">**methods**  
    </span> Starts the METHODS block

-   <span style="color: red">**endmethods**  
    </span> Terminates the METHODS block.

-   <span style="color: red">**method \(mode\)**  
    </span> Sets the algorithm to be executed. The relevant choice for
    \(mode\) is VERTEX FCIMC to run an FCIQMC calculation. Alternative
    choices are DETERM-PROJ to run a deterministic calculation and
    SPECTRAL-LANCZOS to calculate a spectrum using the lanczos
    algorithm.
