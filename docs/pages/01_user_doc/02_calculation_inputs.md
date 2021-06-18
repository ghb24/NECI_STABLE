---
title: Calculation inputs
---

## Calculation inputs

The NECI executable takes one input argument, which is the name of an
input file containing the instructions for carrying out the calculation.
The input file is organized in blocks, with each block being started and
terminated by a dedicated keyword. Each block can contain a number of
keywords to specify options. Here, a list of the blocks and their
respective keywords is given.

The first line of the input is always `title`, the last line is always
`end`.

Some keywords are mandatory, those are marked in **<span style="color: red">red</span>** and are given at the beginning of the
description of each paragraph. Then come recommended options, marked in
**<span style="color: blue">blue</span>**, followed by further options
given in black.

Keywords which are purely for debugging purposes and only interesting
for developers are markes as **<span style="color: green">green</span>**.

### SYSTEM Block

The SYSTEM block specifies the properties of the physical system that is
considered. The block starts with the `system` keyword and ends with the
`endsys` keyword.

-   **<span style="color: red">system</span>**<br>
    Starts the SYSTEM block. Has one mandatory additional
    argument to specify the type of the system. The options are

    -   **read**<br>
        Read in the integrals from a FCIDUMP file, used for ab-initio
        calculations.

    -   **FCIDUMP-name**<br>
        Specify the name of the FCIDUMP file. Defaults to FCIDUMP.

    -   **hubbard**<br>
        Uses the Hubbard model Hamiltonian.

        -   **k-space**<br>
            In the momentum-space basis (beneficial for low \(U/t\))

        -   **real-space**<br>
            In the real-space basis (beneficial for large \(U/t\))

    -   **ueg**<br>
        Uses the Hamiltonian of the uniform electron gas in a box.

-   **<span style="color: red">endsys</span>**<br>
 Terminates the SYSTEM block.

-   **<span style="color: red">electrons \(n\), nel \(n\)</span>**<br>
 Sets the number of electrons to \(n\)

-   **spin-restrict [\(m\)]**<br>
    Sets the total \(S_z\) quantum number to \(\frac{m}{2}\). The
    argument \(m\) is optional and defaults to 0.

-   **hphf \(s\)**<br>
    Uses a basis of (anti-)symmetric combinations of Slater determinants
    with respect to global spin-flip. \(s=0\) indicates anti-symmetric
    combinations, \(s=1\) symmetric combinations. This is useful to
    exclude unwanted spin configurations. For example, no triplet states
    can occur for `hphf 0`.

-   **guga \(S\)**<br>
    Activates the spin-adapted implementation of FCIQMC using CSFs (more
    precisely Gel’fand-Tsetlin states) and the graphical Unitary Group
    Approach (GUGA). The total spin \(S\) must be specified in multiples
    of \(\hbar/2\). So \(S = 0\) is a singlet, \(S = 1\) a doublet
    \(S = 2\) a Triplet and so on. This keyword MUST be combined with
    certain `nonuniformrandexcits` input as described below! Also it is
    not allowed to use this keyword with the `spin-restrict` or `HPHF`
    keyword. And the number of electrons must allow the chosen spin
    (even number for Singlet,Triplet,etc. and odd for
    Doublet,Quartet,...). Additionally the choice of a reference
    determinant via the `definedet` keyword of the `CALC` input block
    (described below) is necessary in most cases, as the automatic
    reference state setup is not always working in a CSF-based
    implementation.

-   **sym \(k_x\) \(k_y\) \(k_z\) \(s\)**<br>
    Specifies the symmetry of the target state. The first three
    arguments set the momentum \((k_x,k_y,k_z)\) and are only used for
    Hubbard and ueg-type systems, the last argument \(s\) specifies the
    irrep within \(d_{2h}\) and is only used for ab-initio systems.

-   **lztot**<br>
    Set the total \(L_s\) quantum number. Has one mandatory additional
    argument, which is the value of \(L_s\).

-   **useBrillouinTheorem**<br>
    Assume that single excitations have zero matrix elements with the
    reference. By default, this is determined automatically.

-   **noBrillouinTheorem**<br>
    Always assume that single excitations have nozero matrix elements
    with the reference.

-   **umatEpsilon \(\epsilon\)**<br>
    Defines a threshold value \(\epsilon\) below which matrix elements
    of the Hamiltonian are rounded to 0. Defaults to \(10^{-8}\).

-   **diagonaltmat**<br>
    Assume the kinetic operator is diagonal in the given basis set.

-   **noSingExcits**<br>
    Assume there is no coupling between single excitations in the
    Hamiltonian.

-   **rohf**<br>
    Use restricted open-shell integrals.

-   **read\_rofcidump**<br>
    Read the integrals from a ROFCIDUMP file.

-   **spinorbs**<br>
    Uses spin orbitals instead of spatial orbtials for addressing the
    integrals. This can be used if the integrals depend on the spin of
    the orbitals.

-   **molproMimic**<br>
    Use the same orbital ordering as molpro, mimicking the behaviour of
    calling NECI from molpro. First nelec/2 orbitals sorted by diagonal
    elements of the Fock matrix are considered to be occupied, and the
    rest to be virtual. Each sector is then sorted separately by
    symmetry labels.

-   **complexOrbs\_realInts**<br>
    The orbitals are complex, but not the integrals. This reduces the
    symmetry of the 4-index integrals. Only affects `kneci` and `kmneci`
    calculations.

-   **complexWalkers-realInts**<br>
    The integrals and orbitals are real, but the wave function shall be
    complex. Only affects `kneci` and `kmneci` calculations.

-   **system-replicas \(n\)**<br>
    Specifies the number of wave functions that shall be evolved in
    parallel. The argument \(n\) is the number of wave functions
    (replicas). Requires `mneci` or `kmneci`.

#### Excitation generation options

-   **<span style="color: blue">nonUniformRandExcits</span>**<br>
 Use a non-uniform random excitation generator for picking
    the move in the FCIQMC spawn step. This can significantly speed up
    the calculation. Requires an additional argument, that can be chosen
    from the following

    -   **<span style="color: blue">pchb</span>**<br>
 Generates excitations weighted directly with the matrix
        elements using pre-computed alias tables. This excitation
        generator is extremely fast, while maintaining high acceptance
        rates and is generally recommended when memory is not an issue.

    -   **nosymgen**<br>
        Generate all possible excitations, regardless of symmetry. Might
        have a low acceptance rate.

    -   **4ind-weighted**<br>
        Generate exictations weighted by a Cauchy-Schwarz estimate of
        the matrix element. Has very good acceptance rates, but is
        comparably slow. Using the `4ind-weighted-2` or
        `4IND-WEIGHTED-UNBOUND` instead is recommended.

    -   **4ind-weighted-2**<br>
        Generates excitations using the same Cauchy-Schwary estimate as
        `4ind-weighted`, but uses an optimized algorithm to pick
        orbitals of different spin, being faster than the former.

    -   **4ind-weighted-unbound**<br>
        Generates excitations using the same Cauchy-Schwary estimate as
        `4-ind-weighted` and optimizations as `4ind-weighted-2`, but
        uses more accurate estimates, having higher acceptance rates.
        This excitation generator has high acceptance rates at
        negligible memory cost.

    -   **pcpp**<br>
        The pre-computed power-pitzer excitation generator [@Neufeld2019]. Has low
        memory cost and scales only mildly with system size, and can
        thus be used for large systems.

    -   **mol-guga-weighted**<br>
        Excitation generator for molecular systems used in the
        spin-adapted GUGA Approach. The specification of this excitation
        generator when using GUGA in ab initio systems is necessary!

    -   **ueg-guga**<br>
        Excitation generator choice when using GUGA in the UEG or
        k-space/real-space Hubbard model calculations. It is mandatory
        to specify this keyword in this case!

-   **lattice-excitgen**<br>
    Generates uniform excitations using momentum conservation. Requires
    the `kpoints` keyword.

-   **pchb-weighted-singles**<br>
    Use a weighted single excitation generator for the pchb excitation
    generator. By default, singles are created uniformly, the weighted
    generation is much more expensive, but can help if single matrix
    elements are large.

<!-- -->

-   **GAS-SPEC**<br>
    Perform a *Generalized Active Spaces* (GAS) calculation and specify
    the GAS spaces. It is possible to select the actual implementation
    with the `GAS-CI` keyword. It is possible to use *local* or
    *cumulative* constraints on the particle number. Local constraints
    define the minimum and maximum particle number per GAS space.
    Cumulative constraints define cumulative minima and maxima of the
    cumulative particle number. The specification is first `LOCAL` or
    `CUMULATIVE` to define the kind of constraints followed by the
    number of GAS spaces \(n_\text{GAS}\). The next items are
    \(3 \times n_\text{GAS}\) numbers which are the number of spatial
    orbitals and (cumulative) minimum and maximum number of particles
    per GAS space \(n_i, N_i^\text{min}, N_i^\text{max}\). Finally an
    integer array denotes for each spatial orbital to which GAS space it
    belongs. Instead of `1 1 1 1 1` one can write `5*1`. It is
    advantageous to use the line continuation (`+++`) for human-readable
    formatting as table. Two benzenes with single inter-space excitation
    would be e.g. denoted as:

        GAS-SPEC LOCAL 2 +++
                 6  5  7  +++
                 6  5  7  +++
                 1  1  1  1  1  1  2  2  2  2  2  2

    or

        GAS-SPEC LOCAL  2 +++
                 6  5  7  +++
                 6  5  7  +++
                 6*1 6*2

    or

        GAS-SPEC CUMULATIVE 2 +++
                 6  5  7  +++
                 6 12 12 +++
                 6*1 6*2

    In the given example the local and cumulative constraints are
    equivalent, but they are not always!

-   **GAS-CI**<br>
    *Optional keyword.* Specify the actual implementation for GAS. If it
    is ommitted, it will be deduced from `GAS-SPEC`.

    -   **GENERAL-PCHB**<br>
        This is the default and the fastest implementation, if
        sufficient memory is available. The double excitations work
        similar to the FCI precomputed heat bath excitation generator
        but automatically exclude GAS forbidden excitations. If one
        follows the keyword, by `SINGLES`, one can select the singles
        for which there are three possibilities:

        -   **PC-UNIFORM**<br>
            This is the default. It chooses GAS allowed electrons
            uniformly.

        -   **DISCARDING-UNIFORM**<br>
            It chooses electrons uniformly as in FCI and discards.

        -   **ON-FLY-HEAT-BATH**<br>
            It chooses GAS allowed electrons weighted by their matrix
            element.

        An example is

            GAS-CI GENERAL-PCHB +++
                            SINGLES ON-FLY-HEAT-BATH

    -   **DISCARDING**<br>
        Use a Full CI excitation generator and just discard excitations
        which are not contained in the GAS space. Currently PCHB is used
        for Full CI.

    -   **GENERAL**<br>
        Use heat bath on the fly general GAS, which is applicable to any
        GAS specification, but a bit slower than necessary for
        disconnected spaces.

    -   **DISCONNECTED**<br>
        Use the disconnected GAS implementations, which assumes
        disconnected spaces and performs there a bit better than the
        general implementation.

#### Hubbard model and UEG options

-   **lattice \(type\) \(l_x\) \(l_y\) [\(l_Z\)]**<br>
    Defines the basis in terms of a lattice of type \(type\) with extent
    \(l_x
        \times l_y \times l_z\). \(l_x\) and \(l_y\) are mandatory,
    \(l_z\) is optional. \(type\) can be any of `chain`, `star`,
    `square`, `rectangle`, `tilted`, `triangular`, `hexagonal`,
    `kagome`, `ole`.

-   **U \(U\)**<br>
    Sets the Hubbard interaction strengh to \(U\). Defaults to 4.

-   **B \(t\)**<br>
    Sets the Hubbard hopping strengh to \(t\). Defaults to -1.

-   **twisted-bc \(t_1\) [\(t_2\)]**<br>
    Use twisted boundary conditions with a phase of
    \(t_1 \, \frac{2\Pi}{x}\) applied along x-direction, where \(x\) is
    the lattice size. \(t_2\) is optional and the additional phase along
    y-direction, in multiples of \(\frac{2\Pi}{y}\).

-   **open-bc [\(direction\)]**<br>
    Set the boundary condition in \(direction\) to open, i.e. no hopping
    across the cell boundary is possible. \(direction\) is optional and
    can be one of `X`, `Y` or \(\texttt{XY}\), for open boundary
    conditions in x-, y- or both directions. If omitted, both directions
    are given open boundary conditions. Requires a real-space basis.

-   **ueg-offset \(k_x\) \(k_y\) \(k_z\)**<br>
    Offset \((k_x, k_y, k_z)\) for the momentum grid used in for the
    uniform electron gas.

#### Transcorrelation options

-   **molecular-transcorr**<br>
    Enable the usage of a transcorrelated ab-initio Hamiltonian. This
    implies the non-hermiticity of the Hamiltonian as well as the
    presence of 3-body interactions. Requires passing the 3-body
    integrals in either ASCII or HDF5 format in a `TCDUMP` file, or
    `tcdump.h5`, respectively. Enables triple excitation generation.
    When using this option, non-uniform random excitation generator
    become inefficient, so using `nonUniformRandExcits` is discouraged.

-   **ueg-transcorr \(mode\)**<br>
    Enable the usage of a transcorrelated Hamiltonian for the uniform
    electron gas. This implies the non-hermiticity of the Hamiltonian as
    well as the presence of 3-body interactions. \(mode\) can be one of
    `3-body`, `trcorr-excitgen` or `rand-excitgen`.

-   **transcorr [\(J\)]**<br>
    Enable the usage of a transcorrelated Hamiltonian for the real-space
    hubbard mode. This implies the non-hermiticity of the Hamiltonian.
    The optional parameter \(J\) is the transcorrelation parameter and
    defaults to 1.0.

-   **2-body-transcorr [\(J\)]**<br>
    Enable the usage of a transcorrelated Hamiltonian for the momentum
    space hubbard model. This implies the non-hermiticity of the
    Hamiltonian as well as the presence of 3-body interactions. The
    optional argument \(J\) is the correlation parameter and defaults to
    0.25.

-   **exclude-3-body-ex**<br>
    Disables the generation of triple excitations, but still takes into
    account 3-body interactions for all other purposes.

### CALC Block

The CALC block is used to set options concerning the simulation
parameters and modes of FCIQMC. The block starts with the `calc` keyword
and ends with the `endcalc` keyword.

-   **<span style="color: red">calc</span>**<br>
 Starts the CALC block

-   **<span style="color: red">endcalc</span>**<br>
 Terminates the CALC block

-   **<span style="color: blue">time \(t\)</span>**<br>
 Set the maximum time \(t\) in minutes the calculation is
    allowed to run. After \(t\) minutes, the calculation will end.

-   **<span style="color: blue">nmcyc \(n\)</span>**<br>
 Set the maximum number of iterations the calculation is
    allowed to do. After \(n\) iterations, the calculation will end.

-   **eq-cyc \(n\)**<br>
    Set the number of iterations to perform after reaching variable
    shift mode. \(n\) Iterations after entering the variable shift mode,
    the calculation will end. If both `eq-cyc` and `nmcyc` are given,
    the calculation ends when either of the iteration limits is reached.

-   **seed \(s\)**<br>
    Sets the seed of the random number generator to \(s\). This can be
    used to specifically probe for stochastic effects, but is generally
    not required.

-   **averageMcExcits \(x\)**<br>
    Sets the average number of spawning attempts from each walker to
    \(x\).

-   **rdmSamplingIters \(n\)**<br>
    Set the maximum number of iterations used for sampling the RDMs to
    \(n\). After \(n\) iterations of sampling RDMs, the calculation will
    end.

-   **load-balance-blocks [OFF]**<br>
    Distribute the determinants blockwise in a dynamic fashion to
    maintain equal load for all processors. This is enabled by default
    and has one optional argument OFF. If given, the load-balancing is
    disabled.

-   **energy**<br>
    Additionally calculate and print the ground state energy using an
    exact diagonalization technique.

-   **averageMcExcits \(n\)**<br>
    The number of spawns to attempt per walker. Defaults to \(1\) and
    should not be changed without good reason.

-   **adjust-averageMcExcits**<br>
    Dynamically update the number of spawns attempted per walker. Can be
    used if the excitation generator creates a lot of invalid
    excitations, but should be avoided else.

-   **scale-spawns**<br>
    Store the maximum value of \(\frac{H_{ij}}{p_{gen}}\) for each
    determinant and use it to estimate the number of spawns per walker
    to prevent blooms. Useful when this fraction strongly depends on the
    determinant.

#### Population control options

-   **<span style="color: red">totalWalkers \(n\)</span>**<br>
 Sets the targeted number of walkers to \(n\). This means,
    the shift will be varied to keep the walker number constant once it
    reaches \(n\).

-   **<span style="color: blue">diagShift \(S\)</span>**<br>
 Set the initial value of the shift to \(S\). A value of
    \(S<0\) is not recommended, as it will decrease the population from
    the beginning.

-   **<span style="color: blue">shiftDamp \(\zeta\)</span>**<br>
 Set the damping factor used in the shift update scheme to
    \(\zeta\). Defaults to \(10\).

-   **<span style="color: blue">stepsSft \(n\)</span>**<br>
 Sets the number of steps per update cycle of the shift to
    \(n\). Defaults to \(100\).

-   **fixed-n0 \(n_0\)**<br>
    Instead of varying the shift to fix the total number of walkers,
    keep the number of walkers at the reference fixed at \(n_0\).
    Automatically sets `stepsSft 1` and overwrites any `stepssft`
    options given.

-   **targetGrowRate \(grow\) \(walks\)**<br>
    When the number of walkers in the calculation exceeds \(walk\), the
    shift is iteratively adjusted to maintain a fixed grow rate \(grow\)
    until reaching the requested number of total walkers.

-   **jump-shift [OFF]**<br>
    When entering the variable shift mode, the shift will be set to the
    current projected energy. This is enabled by default. There is an
    optional argument OFF that disables this behaviour.

-   **pops-jump-shift**<br>
    Reset the shift when restarting a previous calculation to the
    current projected energy instead of using the shift from the
    previous calculation.

-   **trunc-nopen \(n\)**<br>
    Restrict the Hilbert space of the calculation to those determinants
    with at most \(n\) unpaired electrons.

-   **avGrowthRate [OFF]**<br>
    Average the change in walker number used to calculate the shift.
    This is enabled by default and has one optional argument OFF, which,
    when given, turns the option off.

#### Real walker coefficient options

-   **<span style="color: blue">allRealCoeff</span>**<br>
 Allow determinants to have non-integer population. There is
    a minimal population below which the population of a determinant
    will be rounded stochastically. This defaults to \(1\).

-   **<span style="color: blue">realSpawnCutoff \(x\)</span>**<br>
 Continuous real spawning will be performed, unless the spawn
    has weight less than x. In this case, the weight of the spawning
    will be stochastically rounded up to x or down to zero, such that
    the average weight of the spawning does not change. This is a method
    of removing very low weighted spawnings from the spawned list, which
    require extra memory, processing and communication. A reasonable
    value for x is 0.01.

-   **realCoeffbyExcitLevel \(n\)**<br>
    Allow all determinants up to an excitation level of \(n\) to have
    non-integer population.

-   **setOccupiedThresh \(x\)**<br>
    Set the value for the minimum walker weight in the main walker list.
    If, after all annihilation has been performed, any determinants have
    a total weight of less than x, then the weight will be
    stochastically rounded up to x or down to zero such that the average
    weight is unchanged. Defaults to \(1\), which should only be changed
    with good reason.

-   **energy-scaled-walkers [\(mode\) \(\alpha\) \(\beta\)]**<br>
    Scales the occupied threshold with the energy of a determinant. Has
    three optional arguments and requires `allRealCoeff`. The argument
    \(mode\) can be one of EXPONENTIAL, POWER, EXP-BOUND or NEGATIVE and
    defaults to POWER. Both \(\alpha\) and \(\beta\) default to 1 and
    `realspawncutoff` \(\beta\) is implied.

#### Time-step options

-   **<span style="color: red">tau \(\tau\)</span> [SEARCH]**<br>
    Sets the timestep per iteration to \(\tau\). Has one optional
    argument SEARCH. If given, the time-step will be iteratively updated
    to keep the calculation stable.

-   **<span style="color: blue">hist-tau-search [\(c\) \(nbins\)
    \(bound\)]</span>**<br>
 Update the time-step based on histogramming of the ratio
    \(\frac{H_{ij}}{p(i|j)}\). Not compatible with the tau \(\tau\)
    SEARCH option. The three arguments \(c\), \(nbins\) and \(bound\)
    are optional. \(0<c<1\) is the fraction of the histogram used for
    determining the new timestep, \(nbins\) the number of bins in the
    histogram and \(bound\) is the maximum value of
    \(\frac{H_{ij}}{p(i|j)}\) to be stored.<br>
    For spin-adapted GUGA calculations this option is *highly*
    recommended! Otherwise the time-step can become quite small in these
    simulations.

-   **<span style="color: blue">max-tau \(\tau_\text{max}\)</span>**<br>
 Sets the maximal value of the time-step to
    \(\tau_\text{max}\). Defaults to \(1\).

-   **min-tau [\(\tau_\text{min}\)]**<br>
    Sets the minimal value of the time-step to \(\tau_\text{min}\) and
    enables the iterative update of the time-step. Defaults to
    \(10^{-7}\). The argument \(\tau_\text{min}\) is optional.

-   **keepTauFixed**<br>
    Do never update \(\tau\) and the related parameter
    \(p_\text{singles}\), \(p_\text{doubles}\) or \(p_\text{parallel}\).

-   **truncate-spawns [\(n\) UNOCC]**<br>
    Truncate spawns which are larger than a threshold value \(n\). Both
    arguments are optional, \(n\) defaults to \(3\). If UNOCC is given
    the truncation is restricted to spawns onto unoccupied. Useful in
    combination with hist-tau-search.

-   **maxWalkerBloom \(n\)**<br>
    The time step is scaled such that at most \(n\) walkers are spawned
    in a single attempt, with the scaling being guessed from previous
    spawning attempts.

#### Wave function initialization options

-   **<span style="color: blue">walkContGrow</span>**<br>
 When reading in a wave function from a file, do not set the
    shift or enter variable shift mode.

-   **<span style="color: blue">defineDet \(det\)</span>**<br>
 Sets the reference determinant of the calculation to
    \(det\). If no other initialisation is specified, this will also be
    the initial wave function. The format can either be a
    comma-separated list of spin orbitals, a range of spin orbitals
    (like 12-24) or a combination of both.<br>
    For spin-adapted calculations using the `GUGA` keyword defining a
    starting reference CSF manually is *highly* encouraged, as the
    automatic way often fails. It works in similar ways as for SDs,
    however odd numbered singly occupied orbitals indicate a positive
    spin coupled orbital in GUGA CSFs and even numbered singly occupied
    orbitals negatively spin coupled. So one needs to be careful to not
    define a unphysical CSF with an negative intermediate total spin
    \(S_i < 0\). E.g. a CSF like:<br>
    `definedet 1 2 4 5`<br>
    would cause the calculation to crash, as the first singly occupied
    orbital (4) would cause the total spin to be negative \(S_i < 0\).
    When defining the starting CSF one also needs to ensure that the
    defined CSFs satisfies the total number of electrons and total \(S\)
    defined in the `System` block of the input with the keywords `nel`
    and `guga`. As an example a triplet (`guga 2`) CSF with 4 electrons
    (`nel 4`) would be<br>
    `definedet 1 2 3 5`<br>
    with the first spatial orbital doubly occupied and 2 open-shell
    orbitals (positively spin coupled, hence odd numbers).

-   **readPops**<br>
    Read in the wave function from a file and use the read-in wave
    function for initialisation. In addition to the wave function, also
    the time-step and the shift are read in from the file. This starts
    the calculation in variable shift mode, maintaining a constant
    walker number, unless `walkContGrow` is given.

-   **readpops-changeref**<br>
    Allow the reference determinant to be updated after reading in a
    wave function.

-   **startSinglePart [\(n\)]**<br>
    Initialise the wave function with \(n\) walkers on the reference
    only unless specified differently. The argument \(n\) is optional
    and defaults to \(1\).

-   **proje-changeRef [\(frac\) \(min\)]**<br>
    Allow the reference to change if a determinant obtains \(frac\)
    times the population of the current reference and the latter has a
    population of at least \(min\). Both arguments are optional and
    default to \(1.2\) and \(50\), respectively. This is enabled by
    default.

-   **no-changeref**<br>
    Never change the reference.

#### Initiator options

-   **<span style="color: blue">truncInitiator</span>**<br>
 Use the initiator method [@Cleland2010].

-   **<span style="color: blue">addToInitiator \(x\)</span>**<br>
 Sets the initiator threshold to \(x\), so any determinant
    with more than \(x\) walkers will be an initiator.

-   **senior-initiators [\(age\)]**<br>
    Makes any determinant that has a half-time of at least \(age\)
    iterations an initiator. \(age\) is optional and defaults to \(1\).

-   **superInitiator [\(n\)]**<br>
    Create a list of \(n\) superinitiators, from which all connected
    determinants are set to be initiators. The superinitiators are
    chosen according to population. \(n\) is optional and defaults to
    \(1\).

-   **coherent-superInitiators [\(mode\)]**<br>
    Apply a restriction on the sign-coherence between a determinant and
    any connected superinitiator to determine whether it becomes an
    initiator due to connection. The optional argument \(mode\) can be
    chosen from STRICT, WEAK, XI, AV and OFF. The default is WEAK and is
    enabled by default if the `superInitiators` keyword is given.

-   **dynamic-superInitiators \(n\)**<br>
    Updates the list of superinitiators every \(n\) steps. This is
    enabled by default with \(n=100\) if the `superInitiators` keyword
    is given. A value of 0 indicates no update. This implies the
    `dynamic-core n` option (with the same \(n\)) unless specified
    otherwise.

-   **allow-signed-spawns \(mode\)**<br>
    Never abort spawns with a given sign, regardless of initiators.
    \(mode\) can be either POS or NEG, indicating the sign to keep.

-   **initiator-space**<br>
    Define all determinants within a given initiator space as
    initiators. The space is specified through one of the following
    keywords

    -   **doubles-initiator**<br>
        Use the reference determinant and all single and double
        excitations from it to form the initiator space.

    -   **cas-initiator cas1 cas2**<br>
        Use a CAS space to form the initiator space. The parameter cas1
        specifies the number of electrons in the cas space and cas2
        specifies the number of virtual spin orbitals (the cas2 highest
        energy orbitals will be virtuals).

    -   **ras-initiator ras1 ras2 ras3 ras4 ras5**<br>
        Use a RAS space to form the initiator space. Suppoose the list
        of spatial orbitals are split into three sets, RAS1, RAS2 and
        RAS 3, ordered by their energy. ras1, ras2 and ras3 then specify
        the number of spatial orbitals in RAS1, RAS2 and RAS3. ras4
        specifies the minimum number of electrons in RAS1 orbitals. ras5
        specifies the maximum number of electrons in RAS3 orbitals.
        These together define the RAS space used to form the initiator
        space.

    -   **optimised-initiator**<br>
        Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
        230201). One also needs to use either
        optimised-initiator-cutoff-amp or optimised-initiator-cutoff-num
        with this option.

    -   **optimised-initiator-cutoff-amp \(x1\), \(x2\), \(x3\)...**<br>
        Perform the optimised initiator option, and in iteration \(i\),
        choose which determinants to keep by choosing all determinants
        with an amplitude greater than \(xi\) in the ground state of the
        space (see PRL 109, 230201). The number of iterations is
        determined by the number of parameters provided.

    -   **optimised-initiator-cutoff-num \(n1\), \(n2\), \(n3\)...**<br>
        Perform the optimised initiator option, and in iteration \(i\),
        choose which determinants to keep by choosing the ni most
        significant determinants in the ground state of the space (see
        PRL 109, 230201). The number of iterations is determined by the
        number of parameters provided.

    -   **fci-initiator**<br>
        Use all determinants to form the initiator space. A fully
        deterministic projection is therefore performed with this
        option.

    -   **pops-initiator \(n\)**<br>
        When starting from a POPSFILE, this option will use the \(n\)
        most populated determinants from the popsfile to form the
        initiator space.

    -   **read-initiator**<br>
        Use the determinants in the INITIATORSPACE file to form the
        initiator space. A INITIATORSPACE file can be created by using
        the write-initiator option in the LOGGING block.

#### Adaptive shift options

-   **auto-adpative-shift [\(t\) \(\alpha\) \(c\)]**<br>
    Scale the shift per determinant based on the acceptance rate on a
    determinant. Has three optional arguments. The first is the
    threshold value \(t\) which is the minimal number of spawning
    attempts from a determinant over the full calculation required
    before the shift is scaled, with a default of \(10\). The second is
    the scaling exponent \(\alpha\) with a default of \(1\) and the last
    is the minimal scaling factor, which uses
    \(\frac{1}{\text{HF conn.}}\) as default.

-   **linear-adaptive-shift [\(\sigma\) \(f_1\) \(f_2\)]**<br>
    Scale the shift per determinant linearly with the population of a
    determinant. All arguments are optional and define the function used
    for scaling. \(\sigma\) gives the minimal walker number required to
    have a shift and defaults to \(1\), \(f_1\) the shift fraction to be
    applied at \(\sigma\) with a default of \(0\) and \(f_2\) is the
    shift fraction to be applied at the initiator threshold, defaults to
    \(1\). Every initiator is applied the full shift.

-   **exp-adaptive-shift [\(\alpha\)]**<br>
    Scales the shift expoentially with the population of a determinant.
    The optional argument \(\alpha\) is the exponent of scaling, the
    default is \(2\).

-   **core-adaptive-shift**<br>
    By default, determinants in the corespace are always applied the
    full shift. Using this option also scales the shift in the
    corespace.

-   **aas-matele2**<br>
    Uses the matrix elements for determining the scaling factor in the
    `auto-adaptive-shift`. The recommended option to scale the shift.

#### Multi-replica options

-   **multiple-initial-refs**<br>
    Define a reference determinant for each replica. The following \(n\)
    lines give the reference determinants as comma-separated lists of
    orbitals, where \(n\) is the number of replicas.

-   **orthogonalise-replicas**<br>
    Orthogonalise the replicas after each iteration using Gram Schmidt
    orthogonalisation. This will converge each replica to another state
    in a set of orthogonal eigenstates. Can be used for excited state
    search.

-   **orthogonalise-replicas-symmetric**<br>
    Use the symmetric Löwdin orthonaliser instead of Gram Schmidt for
    orthogonalising the replicas.

-   **replica-single-det-start**<br>
    Starts each replica from a different excited determinant.

#### Semi-stochastic options

-   **<span style="color: blue">semi-stochastic</span>**<br>
 Turn on the semi-stochastic adaptation.

-   **<span style="color: blue">pops-core \(n\)</span>**<br>
 This option will use the \(n\) most populated determinants
    to form the core space.

-   **doubles-core**<br>
    Use the reference determinant and all single and double excitations
    from it to form the core space.

-   **cas-core cas1 cas2**<br>
    Use a CAS space to form the core space. The parameter cas1 specifies
    the number of electrons in the cas space and cas2 specifies the
    number of virtual spin orbitals (the cas2 highest energy orbitals
    will be virtuals).

-   **ras-core ras1 ras2 ras3 ras4 ras5**<br>
    Use a RAS space to form the core space. Suppoose the list of spatial
    orbitals are split into three sets, RAS1, RAS2 and RAS 3, ordered by
    their energy. ras1, ras2 and ras3 then specify the number of spatial
    orbitals in RAS1, RAS2 and RAS3. ras4 specifies the minimum number
    of electrons in RAS1 orbitals. ras5 specifies the maximum number of
    electrons in RAS3 orbitals. These together define the RAS space used
    to form the core space.

-   **optimised-core**<br>
    Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
    230201). One also needs to use either optimised-core-cutoff-amp or
    optimised-core-cutoff-num with this option.

-   **optimised-core-cutoff-amp \(x1\), \(x2\), \(x3\)...**<br>
    Perform the optimised core option, and in iteration \(i\), choose
    which determinants to keep by choosing all determinants with an
    amplitude greater than \(xi\) in the ground state of the space (see
    PRL 109, 230201). The number of iterations is determined by the
    number of parameters provided.

-   **optimised-core-cutoff-num \(n1\), \(n2\), \(n3\)...**<br>
    Perform the optimised core option, and in iteration \(i\), choose
    which determinants to keep by choosing the ni most significant
    determinants in the ground state of the space (see PRL 109, 230201).
    The number of iterations is determined by the number of parameters
    provided.

-   **fci-core**<br>
    Use all determinants to form the core space. A fully deterministic
    projection is therefore performed with this option.

-   **read-core**<br>
    Use the determinants in the CORESPACE file to form the core space. A
    CORESPACE file can be created by using the write-core option in the
    LOGGING block.

-   **dynamic-core \(n\)**<br>
    Update the core space every \(n\) iterations, where \(n\) is
    optional and defaults to \(400\). This is enabled by default if the
    `superinitiators` option is given.

#### Trial wave function options

-   **<span style="color: blue">trial-wavefunction [\(n\)]</span>**<br>
 Use a trial wave function to obtain an estimate for the
    energy, as described in
    <a href="#sec:trial" data-reference-type="ref" data-reference="sec:trial">1.4</a>.
    The argument \(n\) is optional, when given, the trial wave function
    will be initialised \(n\) iterations after the variable shift mode
    started, else, at the start of the calculation. The trial wave
    function is defined through one of the following keywords

    -   **<span style="color: blue">pops-trial \(n\)</span>**<br>
 When starting from a POPSFILE, this option will use the
        \(n\) most populated determinants from the popsfile to form the
        trial space.

    -   **doubles-trial**<br>
        Use the reference determinant and all single and double
        excitations from it to form the trial space.

    -   **cas-trial cas1 cas2**<br>
        Use a CAS space to form the trial space. The parameter cas1
        specifies the number of electrons in the cas space and cas2
        specifies the number of virtual spin orbitals (the cas2 highest
        energy orbitals will be virtuals).

    -   **ras-trial ras1 ras2 ras3 ras4 ras5**<br>
        Use a RAS space to form the trial space. Suppoose the list of
        spatial orbitals are split into three sets, RAS1, RAS2 and RAS
        3, ordered by their energy. ras1, ras2 and ras3 then specify the
        number of spatial orbitals in RAS1, RAS2 and RAS3. ras4
        specifies the minimum number of electrons in RAS1 orbitals. ras5
        specifies the maximum number of electrons in RAS3 orbitals.
        These together define the RAS space used to form the trial
        space.

    -   **optimised-trial**<br>
        Use the iterative approach of Petruzielo *et al.* (see PRL, 109,
        230201). One also needs to use either optimised-trial-cutoff-amp
        or optimised-trial-cutoff-num with this option.

    -   **optimised-trial-cutoff-amp \(x1\), \(x2\), \(x3\)...**<br>
        Perform the optimised trial option, and in iteration \(i\),
        choose which determinants to keep by choosing all determinants
        with an amplitude greater than \(xi\) in the ground state of the
        space (see PRL 109, 230201). The number of iterations is
        determined by the number of parameters provided.

    -   **optimised-trial-cutoff-num \(n1\), \(n2\), \(n3\)...**<br>
        Perform the optimised trial option, and in iteration \(i\),
        choose which determinants to keep by choosing the ni most
        significant determinants in the ground state of the space (see
        PRL 109, 230201). The number of iterations is determined by the
        number of parameters provided.

    -   **fci-trial**<br>
        Use all determinants to form the trial space. A fully
        deterministic projection is therefore performed with this
        option.

#### Memory options

-   **memoryFacPart \(x\)**<br>
    Sets the factor between the allocated space for the wave function
    and the required memory for the specified number of walkers to
    \(x\). Defaults to \(10\).

-   **memoryFacSpawn \(x\)**<br>
    Sets the factor between the allocated space for new spawns and the
    estimate of required memory for the spawns of the specified number
    of walkers on a single processor to \(x\). The memory required for
    spawns increases, the more processors are used, so when running with
    few walkers on relatively many processors, a large factor might be
    needed. Defaults to \(3\).

-   **prone-walkers**<br>
    Instead of terminating when running out of memory, randomly delete
    determinants with low population and few spawns.

-   **store-dets**<br>
    Employ extra memory to store additional information on the
    determinants that had to be computed on the fly else. Trades in
    memory for faster iterations.

#### Reduced density matrix (RDM) options

-   **rdmSamplingIters \(n\)**<br>
    Set the number of iterations for sampling the RDMs to \(n\). After
    \(n\) iterations of sampling, the calculation ends.

-   **inits-rdm**<br>
    Only take into account initiators when calculating RDMs. By default,
    only restricts to initiators for the right vector used in RDM
    calculation. This makes the RDMs non-variational, and the resulting
    energy is the projected energy on the initiator space.

-   **strict-inits-rdm**<br>
    Require both sides of the inits-rdm to be initiators.

-   **no-lagrangian-rdms**<br>
    This option disables the correction used for RDM calculation for the
    adaptive shift. Use this only for debugging purposes, as the
    resulting RDMs are flawed.

#### METHODS Block

The METHODS block is a subblock of CALC, i.e. it is specified inside the
CALC block. It sets the main algorithm to be used in the calculation.
The subblock is started with the `methods` keyword and terminated with
the `endmethods` keyword.

-   **<span style="color: red">methods</span>**<br>
 Starts the METHODS block

-   **<span style="color: red">endmethods</span>**<br>
 Terminates the METHODS block.

-   **<span style="color: red">method \(mode\)</span>**<br>
 Sets the algorithm to be executed. The relevant choice for
    \(mode\) is VERTEX FCIMC to run an FCIQMC calculation. Alternative
    choices are DETERM-PROJ to run a deterministic calculation and
    SPECTRAL-LANCZOS to calculate a spectrum using the lanczos
    algorithm.

### INTEGRAL Block

The INTEGRAL block can be used to freeze orbitals and set properties of
the integrals. The block is started with the `integral` keyword and
terminated with the `endint` keyword.

-   **integral**<br>
    Starts the INTEGRAL block.

-   **endint**<br>
    Terminates the INTEGRAL block.

-   **freeze \(n\) \(m\)**<br>
    Freeze \(n\) core and \(m\) virtual orbitals which are not to be
    considered active in this calculation. The orbitals are selected
    according to orbital energy, the \(n\) lowest and \(m\) highest
    orbitals in energy are frozen.

-   **freezeInner \(n\) \(m\)**<br>
    Freeze \(n\) core and \(m\) virtual orbitals which are not to be
    considered active in this calculation. The orbitals are selected
    according to orbital energy, the \(n\) highest and \(m\) lowest
    orbitals in energy are frozen.

-   **partiallyFreeze \(n_\text{orb}\) \(n_\text{holes}\)**<br>
    Freeze \(n_\text{orb}\) core orbitals partially. This means at most
    \(n_\text{holes}\) holes are now allowed in these orbitals.

-   **partiallyFreezeVirt \(n_\text{orb}\) \(n_\text{els}\)**<br>
    Freeze \(n_\text{orb}\) virtual orbitals partially. This means at
    most \(n_\text{els}\) electrons are now allowed in these orbitals.

-   **hdf5-integrals**<br>
    Read the 3-body integrals for a transcorrelated ab-initio
    Hamiltonian from an HDF5 file. Ignored when the
    `molecular-transcorrelated` keyword is not given. Requires compiling
    with HDF5.

-   **sparse-lmat**<br>
    Store the 3-body integrals in a sparse format to save memory.
    Initialisation and iterations might be slower. Requires
    `hdf5-integrals`.

### Development keywords

-   **<span style="color: green">UNIT-TEST-PGEN</span>**<br>
    Test the pgens for the \(n(\text{most populated})\) configurations
    of an existing pops-file. Perform the tests \(n(\text{iterations})\)
    times for each configuration. The order of arguments is
    \(n(\text{most populated}), n(\text{iterations})\).

### KP-FCIQMC Block

This block enables the Krylov-projected FCIQMC (KPFCIQMC) method [@Blunt2015]
which is fully implemented in NECI. It requires `dneci` or `mneci` to be
run. When specifying the KP-FCIQMC block, the METHODS block should be
omitted. This block is started with the `kp-fciqmc` keyword and
terminated with the `end-kp-fciqmc` keyword.

-   **kp-fciqmc**<br>
    Starts the KP-FCIQMC block

-   **end-kp-fciqmc**<br>
    Terminates the KP-FCIQMC block.

-   **num-krylov-vecs \(N\)**<br>
    \(N\) specifies the total number of Krylov vectors to sample.

-   **num-iters-between-vecs \(N\)**<br>
    \(N\) specifies the (constant) number of iterations between each
    Krylov vector sampled. The first Krylov vector is always the
    starting wave function.

-   **num-iters-between-vecs-vary \(i_{12}\), \(i_{23}\),
    \(i_{34}\)...**<br>
    \(i_{n,n+1}\) specifies the number of iterations between the nth and
    (n+1)th Krylov vectors. The number of parameters input should be the
    number of Krylov vectors asked for minus one. The first Krylov
    vector is always the starting wave function.

-   **num-repeats-per-init-config \(N\)**<br>
    \(N\) specifies the number repeats to perform for each initial
    configuration, i.e. the number of repeats of the whole evolution,
    from the first sampled Krylov vector to the last. The projected
    Hamiltonian and overlap matrix estimates will be output for each
    repeat, and the averaged values of these matrices used to compute
    the final results.

-   **averagemcexcits-hamil \(N\)**<br>
    When calculating the projected Hamiltonian estimate, an FCIQMC-like
    spawning is used, rather than calculating the elements exactly,
    which would be too computationally expensive. Here, \(N\) specifies
    the number of spawnings to perform from each walker from each Krylov
    vector when calculating this estimate. Thus, increasing \(N\) should
    improve the quality of the Hamiltonian estimate.

-   **finite-temperature**<br>
    If this option is included then a finite-temperature calculation is
    performed. This involves starting from several different random
    configurations, whereby walkers are distributed on random
    determinants. The number of initial configurations should be
    specified with the num-init-configs option.

-   **num-init-configs \(N\)**<br>
    \(N\) specifies the number of initial configurations to perform the
    sampling over. An entire FCIQMC calculation will be performed, and
    an entire subspace generated, for each of these configurations. This
    option should be used with the finite-temperature option, but is not
    necessary for spectral calculations where one always starts from the
    same initial vector.

-   **memory-factor \(x\)**<br>
    This option is used to specify the size of the array allocated for
    storing the Krylov vectors. The number of slots allocated to store
    unique determinants in the array holding all Krylov vectors will be
    equal to \(ABx\), where here \(A\) is the length of the main walker
    list, \(B\) is the number of Krylov vectors, and \(x\) is the value
    input with this option.

-   **num-walker-per-site-init \(x\)**<br>
    For finite-temperature jobs, \(x\) specifies the number of walkers
    to place on a determinant when it is chosen to be occupied.

-   **exact-hamil**<br>
    If this option is specified then the projected Hamiltonian will be
    calculated exactly for each set of Krylov vectors sampled, rather
    than randomly sampling the elements via an FCIQMC-like spawning
    dynamic.

-   **fully-stochastic-hamil**<br>
    If this option is specified then the projected Hamiltonian will be
    estimated without using the semi-stochastic adaptation. This will
    decrease the quality of the estimate, but may be useful for
    debugging or analysis of the method.

-   **init-correct-walker-pop**<br>
    For finite-temperature calculations on multiple cores, the initial
    population may not be quite as requested. This is because the
    quickest (and default) method involves generating determinants
    randomly and sending them to the correct processor at the end. It is
    possible in this process that walkers will die in annihilation.
    However, if this option is specified then each processor will throw
    away spawns to other processors, thus allowing the correct total
    number of walkers to be spawned.

-   **init-config-seeds seed1, seed2...**<br>
    If this option is used then, for finite-temperature calculations, at
    the start of each calculation over an initial configuration, the
    random number generator will be re-initialised with the
    corresponding input seed. The number of seeds provided should be
    equal to the number of initial configurations.

-   **all-sym-sectors**<br>
    If this option is specified then the FCIQMC calculation will be run
    in all symmetry sectors simultaneously. This is an option relevant
    for finite-temperature calculations.

-   **scale-population**<br>
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

### REALTIME Block

The REALTIME block enables the calculation of the real-time evolution of
a given state and is only required for real-time calculations. Real-time
evolution strictly requires `kneci` or `kmneci` and using `kmneci` is
strongly recommended. This block is started with the `realtime` keyword
and terminatd with the `endrealtime` keyword.

-   **realtime**<br>
    Starts the realtime block. This automatically enables `readpops`,
    and reads from a file named `<popsfilename>.0`.

-   **endrealtime**<br>
    Terminates the REALTIME block.

-   **single \(i\) \(a\)**<br>
    Applies a single excitation operator exciting from orbital \(i\) to
    orbital \(a\) to the initial state before starting the calculation.

-   **lesser \(i\) \(j\)**<br>
    Calculates the one-particle Green’s function \(<c_i^\dagger(t)
          c_j>\). Requires using one electron less than in the popsfile.

-   **greater \(i\) \(j\)**<br>
    Calculates the one-particle Green’s function \(<c_i(t)
          c_j^\dagger>\). Requires using one electron more than in the
    popsfile.

-   **start-hf**<br>
    Do not read in a popsfile but start from a single determinant.

-   **rotate-time \(\alpha\)**<br>
    Calculates the evolution along a trajectory \(e^{i\alpha}t\) instead
    of pure real-time trajectory.

-   **dynamic-rotation [\(\eta\)]**<br>
    Determine a time-dependent \(\alpha\) on the fly for `rotate-time`.
    \(\eta\) is optional and a damping parameter, defaults to 0.05.
    \(\alpha\) is then chosen such that the walker number remains
    constant.

-   **rotation-threshold \(N\)**<br>
    Grow the population to \(N\) walkers before starting to adjust
    \(\alpha\).

-   **stepsalpha \(n\)**<br>
    The number of timesteps between two updates of the \(\alpha\)
    parameter when using `dynamic-rotation`.

-   **log-trajectory**<br>
    Output the time trajectory to a separate file. This is useful if the
    same calculation shall be reproduced.

-   **read-trajectory**<br>
    Read the time trajectory from disk, using a file created by NECI
    with the `log-trajectory` keyword.

-   **live-trajectory**<br>
    Read the time trajectory from disk, using a file which is currently
    being created by NECI with the `log-trajectory` keyword. Can be used
    to create additional data for the same trajectory while the original
    calculation is still running.

-   **noshift**<br>
    Do not apply a shift during the real-time evolution. Strongly
    recommended.

-   **stabilize-walkers [\(S\)]**<br>
    Use the shift to stabilize the walker number if it drops below 80%
    of the peak value. \(S\) Is optional and is an asymptotic value used
    to fix the shift.

-   **energy-benchmark \(E\)**<br>
    Set the energy origin to \(E\) by applying a global, constant shift
    to the Hamiltonian. Can be chosen arbitrarily, but a reasonable
    selection can greatly help efficiency.

-   **rt-pops**<br>
    A second popsfile is supplied containing a time evolved state
    created with the `realtime` keyword whose time evolution is to be
    continued. In this case, the original popsfile is still required for
    calculating properties, so two popsfiles will be read in.

### LOGGING Block

The LOGGING block specifies the output of the calculation and which
status information of the calculation shall be collected. This block is
started with the `logging` keyword and terminated with the `endlog`
keyword.

-   **<span style="color: blue">logging</span>**<br>
 Starts the LOGGING block.

-   **<span style="color: blue">endlog</span>**<br>
 Terminates the LOGGING block.

-   **<span style="color: blue">hdf5-pops</span>**\
 Sets the format to read and write the wave function to HDF5.
    Requires building with the `ENABLE-HDF5` cmake option.

-   **popsfile \(n\)**<br>
    Save the current wave function on disk at the end of the
    calculation. Can be used to initialize subsequent calculations and
    continue the run. This is enabled by default. \(n\) is optional and,
    when given, specifies that every \(n\) iteration, the wave function
    shall be saved. Setting \(n=-1\) disables this option.

-   **popsFileTimer \(n\)**<br>
    Write out a the wave function to disk every \(n\) minutes, each time
    overwriting the last output.

-   **hdf5-pops-write**<br>
    Sets the format to write the wave function to HDF5. Requires
    building with the `ENABLE-HDF5` cmake option.

-   **hdf5-pops-read**<br>
    Sets the format to read the wave function to HDF5. Requires building
    with the `ENABLE-HDF5` cmake option.

-   **highlyPopWrite \(n\)**<br>
    Print out the \(n\) most populated determinants at the end of the
    calculation. Is enabled by default with \(n=15\).

-   **replicas-popwrite [\(n\)]**<br>
    Have each replica print out its own highly populated determinants in
    a separate block instead of one collective output of the on average
    most important determinants. Optionally, \(n\) is the numbers of
    determinants to print, this is the same \(n\) as in
    `highlypopwrite`.

-   **inits-exlvl-write \(n\)**<br>
    Sets the excitation level up to which the number of initiators is
    logged to \(n\). Defaults to \(n=8\).

-   **binarypops**<br>
    Sets the format to write the wave function to binary.

-   **nomcoutput**<br>
    Suppress the printing of iteration information to stdout. This data
    is still written to disk.

-   **stepsOutput \(n\)**<br>
    Write the iteration data to stdout/disk every \(n\) iterations.
    Defaults to the number of iterations per shift cycle. Setting
    \(n=0\) disables iteration data output to stdout and uses the shift
    cycle for disk output.

-   **fval-energy-hist**<br>
    Create a histogram of the scaling factor used for the auto-adaptive
    shift over the energy of a determinant. Only has an effect if
    `auto-adaptive-shift` is used.

-   **fval-pops-hist**<br>
    Create a histogram of the scaling factor used for the auto-adaptive
    shift over the population of a determinant. Only has an effect if
    `auto-adaptive-shift` is used.

#### Semi-stochastic output options

-   **write-core**<br>
    When performing a semi-stochastic calculation, adding this option to
    the Logging block will cause the core space determinants to be
    written to a file called CORESPACE. These can then be further read
    in and used in subsequent semi-stochastic options using the
    `read-core` option in the CALC block.

-   **write-most-pop-core-end \(n\)**<br>
    At the end of a calculation, output the \(n\) most populated
    determinants to a file called CORESPACE. This can further be read in
    and used as the core space in subsequent calculations using the
    `read-core` option.

#### RDM output options

These options control how the RDMs are printed. For a description of how
the RDMs are calculated and the content of the files, please see section
<a href="#sec:rdms" data-reference-type="ref" data-reference="sec:rdms">1.7</a>.

-   **calcRdmOnfly \(i\) \(step\) \(start\)**<br>
    Calculate RDMs stochastically over the course of the calculation.
    Starts sampling RDMs after \(start\) iterations, and outputs an
    average every \(step\) iterations. \(i\) indicates whether only
    1-RDMs (1), only 2-RDMs (2) or both are produced.

-   **rdmLinSpace \(start\) \(n\) \(step\)**<br>
    A more user friendly version of `calcrdmonfly` and
    `rdmsamplingiters`, this samples both 1- and 2-RDMs starting at
    iteration \(start\), outputting an average every \(step\) iterations
    \(n\) times, then ending the calculation.

-   **diagFlyOneRdm**<br>
    Diagonalise the 1-RDMs, yielding the occupation numbers of the
    natural orbitals.

-   **printOneRdm**<br>
    Always output the 1-RDMs to a file, regardless of which RDMs are
    calculated. May compute the 1-RDMs from the 2-RDMs.

-   **writeRdmsToRead off**<br>
    The presence of this keyword overrides the default. If the OFF word
    is present, the unnormalised `TwoRDM_POPS_a***` files will
    definitely not be printed, otherwise they definitely will be,
    regardless of the state of the `popsfile/binarypops` keywords.

-   **readRdms**<br>
    This keyword tells the calculation to read in the `TwoRDM_POPS_a***`
    files from a previous calculation. The restarted calc then continues
    to fill these RDMs from the very first iteration regardless of the
    value put with the `calcRdmOnFly` keyword. The calculation will
    crash if one of the `TwoRDM_POPS_a***` files are missing. If the
    `readRdms` keyword is present, but the calc is doing a
    `StartSinglePart` run, the `TwoRDM_POPS_a***` files will be ignored.

-   **noNormRdms**<br>
    This will prevent the final, normalised `TwoRDM_a***` matrices from
    being printed. These files can be quite large, so if the calculation
    is definitely not going to be converged, this keyword may be useful.

-   **writeRdmsEvery \(iter\)**<br>
    This will write the normalised `TwoRDM_a***` matrices every \(iter\)
    iterations while the RDMs are being filled. At the moment, this must
    be a multiple of the frequency with which the energy is calculated.
    The files will be labelled with incrementing values -
    `TwoRDM_a***.1` is the first, and then next `TwoRDM_a***.2` etc.

-   **write-spin-free-rdm**<br>
    Output the spin-free 2-RDMs to disk at the end of the calculation.

-   **printRoDump**<br>
    Output the integrals of the natural orbitals to a file.

-   **print-molcas-rdms**<br>
    It is now possible to calculate stochastic spin-free RDMs with the
    GUGA implementation. This keyword is necessary if one intends to use
    this feature in conjunction with `Molcas` to perform a spin-free
    Stochastic-CASSCF. It produces the three files `DMAT, PAMAT` and
    `PSMAT`, which are read-in by `Molcas`.
