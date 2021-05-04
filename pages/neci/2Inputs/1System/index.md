title: SYSTEM Block 
--- 

### SYSTEM Block

[TOC]

The SYSTEM block specifies the properties of the physical system that is
considered. The block starts with the `system` keyword and ends with the
`endsys` keyword.

-   <span style="color: red">**system**</span>    
    Starts the SYSTEM block. Has one mandatory additional
    argument to specify the type of the system. The options are

    -   **read**  
        Read in the integrals from a FCIDUMP file, used for ab-initio
        calculations.

    -   **FCIDUMP-name**  
        Specify the name of the FCIDUMP file. Defaults to FCIDUMP.

    -   **hubbard**  
        Uses the Hubbard model Hamiltonian.

        -   **k-space**  
            In the momentum-space basis (beneficial for low \(U/t\))

        -   **real-space**  
            In the real-space basis (beneficial for large \(U/t\))

    -   **ueg**  
        Uses the Hamiltonian of the uniform electron gas in a box.

-   <span style="color: red">**endsys**  
    </span> Terminates the SYSTEM block.

-   <span style="color: red">**electrons \(n\), nel \(n\)**  
    </span> Sets the number of electrons to \(n\)

-   **spin-restrict [\(m\)]**  
    Sets the total \(S_z\) quantum number to \(\frac{m}{2}\). The
    argument \(m\) is optional and defaults to 0.

-   **hphf \(s\)**  
    Uses a basis of (anti-)symmetric combinations of Slater determinants
    with respect to global spin-flip. \(s=0\) indicates anti-symmetric
    combinations, \(s=1\) symmetric combinations. This is useful to
    exclude unwanted spin configurations. For example, no triplet states
    can occur for `hphf 0`.

-   **guga \(S\)**  
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

-   **sym \(k_x\) \(k_y\) \(k_z\) \(s\)**  
    Specifies the symmetry of the target state. The first three
    arguments set the momentum \((k_x,k_y,k_z)\) and are only used for
    Hubbard and ueg-type systems, the last argument \(s\) specifies the
    irrep within \(d_{2h}\) and is only used for ab-initio systems.

-   **lztot**  
    Set the total \(L_s\) quantum number. Has one mandatory additional
    argument, which is the value of \(L_s\).

-   **useBrillouinTheorem**  
    Assume that single excitations have zero matrix elements with the
    reference. By default, this is determined automatically.

-   **noBrillouinTheorem**  
    Always assume that single excitations have nozero matrix elements
    with the reference.

-   **umatEpsilon \(\epsilon\)**  
    Defines a threshold value \(\epsilon\) below which matrix elements
    of the Hamiltonian are rounded to 0. Defaults to \(10^{-8}\).

-   **diagonaltmat**  
    Assume the kinetic operator is diagonal in the given basis set.

-   **noSingExcits**  
    Assume there is no coupling between single excitations in the
    Hamiltonian.

-   **rohf**  
    Use restricted open-shell integrals.

-   **read\_rofcidump**  
    Read the integrals from a ROFCIDUMP file.

-   **spinorbs**  
    Uses spin orbitals instead of spatial orbtials for addressing the
    integrals. This can be used if the integrals depend on the spin of
    the orbitals.

-   **molproMimic**  
    Use the same orbital ordering as molpro, mimicking the behaviour of
    calling NECI from molpro. First nelec/2 orbitals sorted by diagonal
    elements of the Fock matrix are considered to be occupied, and the
    rest to be virtual. Each sector is then sorted separately by
    symmetry labels.

-   **complexOrbs\_realInts**  
    The orbitals are complex, but not the integrals. This reduces the
    symmetry of the 4-index integrals. Only affects `kneci` and `kmneci`
    calculations.

-   **complexWalkers-realInts**  
    The integrals and orbitals are real, but the wave function shall be
    complex. Only affects `kneci` and `kmneci` calculations.

-   **system-replicas \(n\)**  
    Specifies the number of wave functions that shall be evolved in
    parallel. The argument \(n\) is the number of wave functions
    (replicas). Requires `mneci` or `kmneci`.

#### Excitation generation options

-   <span style="color: blue">**nonUniformRandExcits**  
    </span> Use a non-uniform random excitation generator for picking
    the move in the FCIQMC spawn step. This can significantly speed up
    the calculation. Requires an additional argument, that can be chosen
    from the following

    -   <span style="color: blue">**pchb**  
        </span> Generates excitations weighted directly with the matrix
        elements using pre-computed alias tables. This excitation
        generator is extremely fast, while maintaining high acceptance
        rates and is generally recommended when memory is not an issue.

    -   **nosymgen**  
        Generate all possible excitations, regardless of symmetry. Might
        have a low acceptance rate.

    -   **4ind-weighted**  
        Generate exictations weighted by a Cauchy-Schwarz estimate of
        the matrix element. Has very good acceptance rates, but is
        comparably slow. Using the `4ind-weighted-2` or
        `4IND-WEIGHTED-UNBOUND` instead is recommended.

    -   **4ind-weighted-2**  
        Generates excitations using the same Cauchy-Schwary estimate as
        `4ind-weighted`, but uses an optimized algorithm to pick
        orbitals of different spin, being faster than the former.

    -   **4ind-weighted-unbound**  
        Generates excitations using the same Cauchy-Schwary estimate as
        `4-ind-weighted` and optimizations as `4ind-weighted-2`, but
        uses more accurate estimates, having higher acceptance rates.
        This excitation generator has high acceptance rates at
        negligible memory cost.

    -   **pcpp**  
        The pre-computed power-pitzer excitation generator [2]. Has low
        memory cost and scales only mildly with system size, and can
        thus be used for large systems.

    -   **mol-guga-weighted**  
        Excitation generator for molecular systems used in the
        spin-adapted GUGA Approach. The specification of this excitation
        generator when using GUGA in ab initio systems is necessary!

    -   **ueg-guga**  
        Excitation generator choice when using GUGA in the UEG or
        k-space/real-space Hubbard model calculations. It is mandatory
        to specify this keyword in this case!

-   **lattice-excitgen**  
    Generates uniform excitations using momentum conservation. Requires
    the `kpoints` keyword.

-   **pchb-weighted-singles**  
    Use a weighted single excitation generator for the pchb excitation
    generator. By default, singles are created uniformly, the weighted
    generation is much more expensive, but can help if single matrix
    elements are large.

<!-- -->

-   **GAS-SPEC**  
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

-   **GAS-CI**  
    *Optional keyword.* Specify the actual implementation for GAS. If it
    is ommitted, it will be deduced from `GAS-SPEC`.

    -   **GENERAL-PCHB**  
        This is the default and the fastest implementation, if
        sufficient memory is available. The double excitations work
        similar to the FCI precomputed heat bath excitation generator
        but automatically exclude GAS forbidden excitations. If one
        follows the keyword, by `SINGLES`, one can select the singles
        for which there are three possibilities:

        -   **PC-UNIFORM**  
            This is the default. It chooses GAS allowed electrons
            uniformly.

        -   **DISCARDING-UNIFORM**  
            It chooses electrons uniformly as in FCI and discards.

        -   **ON-FLY-HEAT-BATH**  
            It chooses GAS allowed electrons weighted by their matrix
            element.

        An example is

            GAS-CI GENERAL-PCHB +++
                            SINGLES ON-FLY-HEAT-BATH

    -   **DISCARDING**  
        Use a Full CI excitation generator and just discard excitations
        which are not contained in the GAS space. Currently PCHB is used
        for Full CI.

    -   **GENERAL**  
        Use heat bath on the fly general GAS, which is applicable to any
        GAS specification, but a bit slower than necessary for
        disconnected spaces.

    -   **DISCONNECTED**  
        Use the disconnected GAS implementations, which assumes
        disconnected spaces and performs there a bit better than the
        general implementation.

#### Hubbard model and UEG options

-   **lattice \(type\) \(l_x\) \(l_y\) [\(l_Z\)]**  
    Defines the basis in terms of a lattice of type \(type\) with extent
    \(l_x
        \times l_y \times l_z\). \(l_x\) and \(l_y\) are mandatory,
    \(l_z\) is optional. \(type\) can be any of `chain`, `star`,
    `square`, `rectangle`, `tilted`, `triangular`, `hexagonal`,
    `kagome`, `ole`.

-   **U \(U\)**  
    Sets the Hubbard interaction strengh to \(U\). Defaults to 4.

-   **B \(t\)**  
    Sets the Hubbard hopping strengh to \(t\). Defaults to -1.

-   **twisted-bc \(t_1\) [\(t_2\)]**  
    Use twisted boundary conditions with a phase of
    \(t_1 \, \frac{2\Pi}{x}\) applied along x-direction, where \(x\) is
    the lattice size. \(t_2\) is optional and the additional phase along
    y-direction, in multiples of \(\frac{2\Pi}{y}\).

-   **open-bc [\(direction\)]**  
    Set the boundary condition in \(direction\) to open, i.e. no hopping
    across the cell boundary is possible. \(direction\) is optional and
    can be one of `X`, `Y` or \(\texttt{XY}\), for open boundary
    conditions in x-, y- or both directions. If omitted, both directions
    are given open boundary conditions. Requires a real-space basis.

-   **ueg-offset \(k_x\) \(k_y\) \(k_z\)**  
    Offset \((k_x, k_y, k_z)\) for the momentum grid used in for the
    uniform electron gas.

#### Transcorrelation options

-   **molecular-transcorr**  
    Enable the usage of a transcorrelated ab-initio Hamiltonian. This
    implies the non-hermiticity of the Hamiltonian as well as the
    presence of 3-body interactions. Requires passing the 3-body
    integrals in either ASCII or HDF5 format in a `TCDUMP` file, or
    `tcdump.h5`, respectively. Enables triple excitation generation.
    When using this option, non-uniform random excitation generator
    become inefficient, so using `nonUniformRandExcits` is discouraged.

-   **ueg-transcorr \(mode\)**  
    Enable the usage of a transcorrelated Hamiltonian for the uniform
    electron gas. This implies the non-hermiticity of the Hamiltonian as
    well as the presence of 3-body interactions. \(mode\) can be one of
    `3-body`, `trcorr-excitgen` or `rand-excitgen`.

-   **transcorr [\(J\)]**  
    Enable the usage of a transcorrelated Hamiltonian for the real-space
    hubbard mode. This implies the non-hermiticity of the Hamiltonian.
    The optional parameter \(J\) is the transcorrelation parameter and
    defaults to 1.0.

-   **2-body-transcorr [\(J\)]**  
    Enable the usage of a transcorrelated Hamiltonian for the momentum
    space hubbard model. This implies the non-hermiticity of the
    Hamiltonian as well as the presence of 3-body interactions. The
    optional argument \(J\) is the correlation parameter and defaults to
    0.25.

-   **exclude-3-body-ex**  
    Disables the generation of triple excitations, but still takes into
    account 3-body interactions for all other purposes.
