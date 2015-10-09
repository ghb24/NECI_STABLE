.. _input_calc:

----
Calc
----

**CALC**
    Start calculation block.  This chooses what calculation to do.

[Calculation options---see below.]

**ENDCALC**
    End the calculation input block.

General options
---------------

**ALLPATHS**
    Choose all determinants (i.e. set NPATHS = -1).

**BETA** [BETA]
   Set :math:`\beta`.

**BETAOVERP** [BETAP] 
   Default= 1.d-4.

   Set :math:`\beta/P`.

**DELTABETA** [DBETA]
   Set :math:`\delta\beta`.  If given a negative value, calculate it exactly.

   .. note::
     What is this used for?

**DETINV** [DETINV]
    Specify the root determinant for which the complete vertex series is
    worked out, using the determinant index obtained from a previous
    calculation.  If **DETINV** is negative, the NPATHS calculations
    are started at this determinant.

**EXCITE** [ICILEVEL] 
   Default 0.

   Excitiation level at which to truncate determinant list.  If ICILEVEL=0
   then all determinants are enumerated.
   This also works for FCIMC calculations.

**EXCITATIONS** [**OLD** **NEW**]
   For generation of up to double excitations use the old (completely
   reliable), or new (faster, but does not work for more than 2-vertex
   level SUMS) routine

   .. note::
     You can now use the **NEW** routines for all methods, right?
     What is the difference between **NEW** and **OLD**?  (If it doesn't say, how else
     can a user make an informed decision as to which to use?)

**EXCITATIONS** [**SINGLES** **DOUBLES**]
   Default is to use all excitations.

   Restricts determinants which are allowed to be connected to the
   reference determinant to be either single or double excitations of
   the reference determinant.
   
   Applies only to the **VERTEX** [**SUM** **STAR**] **NEW** methods.

**HAMILTONIAN** [**STAR**]
    Store the Hamiltonian.  This is defaulted to ON if **ENERGY** is set,
    but can be used without **ENERGY**.

    **STAR** 
        Only the connections between the root determinant and its
        excitations should be included in the Hamiltonian and not
        off-diagonal elements between excited determinants.

**MAXVERTICES** [MAXVERTICES]
    Give the vertex level of the calculation.  Cannot be used in
    conjunction with a **METHODS** block.

**CONSTRUCTNATORBS**
    Calculates the 1 electron reduced density matrix (1-RDM) as a FCIMC 
    calculation progresses.  At the end of the iterations, this matrix
    is diagonalised to get linear combinations of the HF orbitals which
    approximate the natural orbitals.  The occupation numbers (e-values)
    of these are printed in the OUTPUT file.
    This is now a very old option, a much more efficient equivalent has 
    been added under the **ROTATEORBS** option.  See **USECINATORBS** in the
    system file.


**METHOD** [Method option(s)]
    Specify the method for a graph theory calculation.  See Method
    options for the available methods.

    Can only be specified once if used outside of the methods block, 
    in which case the given method is applied to all vertex levels.

**METHODS**
   Begin a methods block.  This allows a different method for each vertex
   level.  Each vertex level can contain **EXCITATIONS**, **VERTICES**,
   **CYCLES** and **CALCVAR** keywords.
   Each **METHOD** line and the options that follow it detail the calculation
   type for the next vertex level, with the first **METHOD** line used for the 
   the second-vertex level, unless over-ridden with the **VERTICES** option.

   The block terminates with **ENDMETHODS**.

   For example::

      METHODS
         METHOD VERTEX SUM NEW
         EXCITATIONS DOUBLES
         METHOD VERTEX STAR POLY
         EXCITATIONS SINGLES
         VERTICES 2
      ENDMETHODS

   sets the first method, at the two-vertex level, to be a complete 2-vertex
   sum of only doubles, and the second method, overriden to be also at
   the two-vertex level, to be a vertex star of singles.

   Similarly::

      METHODS
         METHOD VERTEX SUM NEW
         METHOD VERTEX SUM MC
         [Monte Carlo options]
      ENDMETHODS

   performs a full sum at the two-vertex level and a Monte Carlo
   calculation at the three-vertex level.

**ENDMETHODS**
   Terminate a methods block.

**PATHS** [option] 
    Select the number of determinants taken to be the root of the graph.
    Usually set to 1.  Valid options:

        NPATHS
            Choose the first NPATHS determinants and calculate RHOPII etc.
        **ALL** 
            Choos all determinants (same as ALLPATHS).
        **ACTIVE** 
            Choose only the active space of determinants: the degenerate
            set containing the highest energy electron.
        **ACTIVE** **ORBITALS** nDown nUp   
            Set the active space to be nDown and nUp orbitals respectively
            from the Fermi level
        **ACTIVE** **SETS** nDown nUp
            Set the active space to be nDown and nUp degenerate sets
            respectively from the Fermi level

**RHOEPSILON** [RHOEPSILON]
    Set the minimum significant value of an element in the :math:`rho`
    matrix as a fraction of the maximum value in the :math:`rho` matrix.
    Matrix elements below this threshold are set to be 0.

**STARCONVERGE** [STARCONV]
    Default 1.d-3.

    Set the convergence criteria for whether a roots to the star graph
    is significant. 

**TROTTER**
   Default.

   Perform a Trotter decomposition to evaluate the :math:`rho` matrix elements.

**TIMESTEPS** [I_P]
    Set P, the timesteps into which :math:`e^{-\beta H}` is split.  Automatically
    sets :math:`\beta/P=0` (as required) but returns an error message if **BETAOVERP** 
    is also used.

**WORKOUT** [NDETWORK]
   Sets the number of determinants which are worked out exactly.

   .. note::
     What is this used for?  

**VERTICES**
   Only available inside a methods block.  
   
   By default, each method takes a
   number of vertices corresponding to its index within the methods
   block, the first methods corresponding to the 2-vertex level, the
   second to the 3-vertex level, and so on.  **VERTICES** overrides this,
   and allows the vertex level of each method to be explicitly specified,
   enabling, for example, the 2-vertex level to be split up and the
   contributions from single and double excitations of the reference
   determinant to be handled separately.

Method options
--------------

**VERTEX SUM** [**OLD** **NEW** **HDIAG**] [**SUB2VSTAR**] [**LOGWEIGHT**]
    Calculate the vertex sum approximation.

    **OLD**
        Diagonalise the :math:`\rho` matrix using the original method.

    **NEW**
        Diagonalise the :math:`\rho` matrix using a more modern, more 
        efficient method.  Recommended.

    **HDIAG**
        Diagonalise the Hamiltonian matrix instead of the :math:`rho` matrix
        in order to calculate the weight and energy contribution of each graph.
    
    **SUB2VSTAR**
        Remove paths which were present in the 2-vertex
        star for each graph.  If this is specified for ANY vertex level,
        it applies to all **SUM** and MC vertex levels.  

    **LOGWEIGHT** 
        Form Q as a multiplication of factors from graphs.  This results
        in the quantity :math:`\operatorname{log} w` being used instead
        of :math:`w`, which also translates to the energy expression
        only involving :math:`\tilde{E}` not weights.  Hopefully this
        is size-consistent.

    .. warning::
      **SUB2VSTAR** and **LOGWEIGHT** are experimental options.

**VERTEX** [**MC** **MCMETROPOLIS** **MCDIRECT** **MCMP**] [**HDIAG**]
    Perform a Monte Carlo calculation.

    **MCDIRECT**
        Perform direct stochastic sampling for the graph theory vertex sum
        method, dividing each freshly generated graph by its normalized
        generation probability.  
        
        If **MULTIMCWEIGHT** is specified then
        the sampling generates graphs from all weighted levels using
        the weighting - a single MC calculation is performed.

        If **MULTIMCWEIGHT** is not specified (default), a separate
        MC calculation is performed at each vertex level.  Combined
        statistics are printed.

        .. warning::
          **MULTIMCWEIGHT** is not documented.  Use with great caution.

    **MCMP**
        Perform direct stochastic sampling, as in **MCDIRECT**,
        but for the Moller--Plesset method.

    **MC** or **MCMETROPOLIS**
        Perform Metropolis Monte Carlo.

        This may be performed in a number of ways. The way is
        chosen by the location of the **VERTEX** **MC** command.

        .. warning:: 

            The following options appear in INPUT_DOC but, however, are incredibly
            poorly documented.  In particular:

                * No detail on the arguments the options take (e.g. **BIAS**).
                * Some options documented don't exist (e.g. **SINGLE**, **BIAS**, **MULTI**, **STOCHASTICTIME**).
                * Sufficient tests are not present in the test suite.

            Do not use.

            The "options" are::

                **STOCHASTICTIME** 
                    may also be specified to perform stochastic
                    time simulations with a given **BIAS**

               **SINGLE**
                   MC is performed at a single vertex level using a composite
                   1-vertex graph containing a full sum previously performed.

               **BIAS** 
                   is used to choose whether a step selects a composite
                   (all lower levels) or a normal (this level) graph.  Stochastic
                   time MC is performed. This can only be specified in the
                   **METHODS** section, and only at the last vertex level.
                   Uses **EXCITWEIGHTING** for excitation generation weighting
                   and **IMPORTANCE** for graph generation weighting

               **MULTI**
                   MC is performed at a multiple vertex levels, but still
                   using a composite 1-vertex graph containing a full sum
                   previously performed. MULTI should be specified in all the
                   (contiguous) vertex levels to be included (not composited)
                   in the MC.  **BIAS** is used to choose whether a step
                   selects a composite (all lower levels) or a normal (the
                   **MULTI** levels) graph.  **MULTIMCWEIGHT** is specified
                   for each **MULTI** level, and gives a relative weighting
                   of selecting the vertex level graphs once a non-composite
                   graph is chosen.  Stochastic time MC is performed.
                   This can only be specified in the **METHODS** section.
                   Once **MULTI** has been specified, it must be specified
                   on all subsequent vertex levels in a **METHODS** section.
                   Uses **EXCITWEIGHTING** for excitation generation weighting
                   and **IMPORTANCE** for graph generation weighting

               **FULL** 
                   Does  MC at all levels using BIAS to bias the levels,
                   **EXCITWEIGHTING** for excitation generation, and
                   **IMPORTANCE** to for graph generation weighting.  This is
                   only available *WITHOUT* a **METHODS** section. If **HDIAG**
                   is specified, the H-diagonalizing routine is used, otherwise,
                   the rho-diagonalizer is used.  **HDIAG** is automatically
                   specified for **MCMP**.

**VERTEX** **SUM** **READ**
    Read in from pre-existing MCPATHS file for that vertex level.
    Only really useful in a **METHODS** section.

**VERTEX** **STAR** [**ADDSINGLES** **COUNTEXCITS**] [star method] [**OLD** **NEW** [**H0**] ] 
    Construct a single and double excitation star from all determinants
    connected to the root (ignoring connections between those dets).
    See [StarPaper]_ for more details.

    **ADDSINGLES** 
        Extend the star graph approach.

        Add the single exctitaions which are en-route to each double
        excitation to that double excitation as spokes, and prediagonalize
        the mini-star centred on each double excitation.  For example,
        if the double excitation is (ij->ab), then singles
        (i->a),(i->b),(j->a) and (j->b) are created in a star with
        (ij->ab), the result diagonalized, and the eigenvalues and
        vectors used to create a new spoke of the main star graph.

        Only works with **NEW**.

    **COUNTEXCITS** 
        Run through all the symmetry allowed excitations
        first and count the connected determinants on the star.  Enables the
        memory requirements to be reduced as only connected determinants need
        to be stored. However, the time taken is increased, as it is necessary
        to run through all determinants in the star twice. Especially useful
        for large systems with memory restraints, when density fitting has
        necessarily turned off symmetry. Also useful if a **RHOEPSILON**
        has been set to a large value so that many of the symmetry allowed
        excitations  will be counted as disconnected.

        .. note::
            Useful for periodic calculations?  Does it need just the
            symmetry info or the transition matrix elements as well?

    **OLD** 
        Use a pre-generated list of determinants using the excitation
        routine version specified in **EXCITATIONS** **OLD** or
        **EXCITATIONS** **NEW**.

    **NEW** 
        Generate determinants on the fly without storing them, using
        the **NEW** excitation routine.  Much more memory efficient.

    **NEW H0** 
        Use the zeroth order N-particle Hamiltonian (shifted such that
        :math:`H^0_{ii} = H_{ii}`) rather than the fully interacting
        Hamiltonian to generate the roots of the polynomial.

        .. note::
          And you'd want to use **NEW H0** why exactly?

    The available star methods are:

        **DIAG** 
            Perform a complete diagonalization on the resultant matrix.  This can
            be very slow. However, by specifying **LANCZOS** in the **CALC**
            block, you can do a Lanczos diagonalisation, which scales much
            better. **EIGENVALUES** can also be specify to only evaluate the
            first few eigenvalues.

        **POLY** 
            Use the special properties of the matrix to find the roots of
            the polynomial and uses them to calculate the relevant values.
            This is order :math:`\text{Ngraph}^2`.

            .. note::
                Ngraph==nDets?

        **POLYMAX** 
            Similar to **POLY** but only finds the highest root of the polynomial, so
            is order Ngraph.  It can be used when P is very large (i.e. :math:`\beta`
            is very large, e.g. 40).

        **POLYCONVERGE** 
            Similar to **POLY** but adds i out of N :math:`\lambda_i`
            roots, such that :math:`(N-i) \lambda_i^P < 10^{-3}`, i.e. we
            evaluate enough roots such that a very conservative error
            estimate of the contribution of the remaining roots is
            negligible.

        **POLYCONVERGE2** 
            Similar to **POLYCONVERGE** but requires 
            :math:`w(1..i) (N-i) \lambda_i^P < 10^{-3}`, where
            :math:`w(1..i)` is the cumulative sum of :math:`\lambda_i^P`,
            which should be a better estimate of the convergence.

    The following are experimental star methods:

        **MCSTAR** 
            Use a basic implementation of the spawning algorithm in
            order to sample the star graph stochastically. The sampling uses
            elements of the Hamiltonian matrix rather than the :math:`rho` matrix, 
            so there will be some differences in the converged energy
            compared to a **VERTEX STAR NEW** calculation.
            
            Many of the **FCIMC** options are also available with MCStar,
            and there are also some extra one.

        **NODAL** 
            Prediagonalise a completely connected set of virtuals for each
            set of occupied (i,j) spin-orbitals. The diagonalised
            excitations are then solved as a star graph. Must be used
            with **NEW**.

        **STARSTARS** 
            Use an approximation that the change of eigenvalues and the
            first element of the eigenvectors of the star graph is linear with
            respect to multiplying the diagonal elements by a constant. Once
            this scaling is found, all stars of stars are prediagonalised,
            and reattached to the original graph. This results in N^2 scaling,
            where N is the number of excitations.

        **TRIPLES** 
            Prediagonalise an excited star of triple excitations from each
            double excitation, reattach the eigenvectors, and solves
            the complete star. Currently only available with '**NEW**',
            '**COUNTEXCITS**' and '**DIAG**'.

Experimental methods
^^^^^^^^^^^^^^^^^^^^

**VERTEX** **FCIMC** [**MCDIFFUSION**] [**RESUMFCIMC**] [**SERIAL**]
    Perform Monte Carlo calculations over pure determinant space, which
    is sampled using a series of 'particles' (or 'walkers').

    The walkers are not necessarily unique and must be sorted at every
    iteration.  Each walker has its own excitation generator.

    **MCDIFFUSION** is a completely particle-conserving diffusion
    algorithm and is much more experimental.

    **FCIMC** and **MCDETS** calculations share many of the same options
    (see Walker Monte Carlo options, below).

    **RESUMFCIMC** creates graphs out of connected determinants, and applies
    the H-matrix successively in order to achieve a local spawning algorithm.
    This reduces to the original spawning algorithm when **GRAPHSIZE** is 2 and
    **HAPP** is 1. Uses many of the same options as **FCIMC**.

    **SERIAL** will force NECI to run the serial FCIMC code (which differs
    substantially from the parallel) even if the code was compiled in parallel.

**VERTEX** **CCMC** [**FCI**] [**EXACTCLUSTER**] [**AMPLITUDE**] [**EXACTSPAWN**] [**BUFFER**]
    Perform Monte Carlo calculations over coupled cluster excitation space, which
    is sampled using a series of 'particles' (or 'walkers').

    The walkers are not necessarily unique and must be sorted at every
    iteration.  Each walker has its own excitation generator.
    **DIRECTANNIHILATION** (in CALC) and **NONUNIFORMRANDEXCITS** (in the SYSTEM section)
    must also be specified.
    
    If **FCI** is specified, then the code runs an equivalent of the **VERTEX** **FCIMC**
    for testing (by only allowing clusters of up to a single excitor, and using (1+T)|D_0>
    as the wavefunction rather than exp(T)|D_0>

    **EXACTCLUSTER** is an exponentially scaling (with number of walkers) algorithm for testing
    the stochastic sampling, and explicitly attempts spawning from all clusters in the space.

    **AMPLITUDE** will enumerate the whole of the allowed space, and assign a floating-point
    amplitude to each excitor.  These amplitudes are stochastically sampled (**INITWALKERS**
    times per MC cycle), and used to propagate the CCMC.

    **EXACTSPAWN** causes spawning to be done exactly - i.e. all allowed connected determinants
    from a given cluster are spawned to.

    **BUFFER** will accumulate all collapsed cluster selections in a buffer and then do spawnings from that.
      When using **EXACTCLUSTER** this is much more efficient.

    Extremely experimental.


**VERTEX** **GRAPHMORPH** [**HDIAG**]
    Set up an initial graph and systematically improve it, by applying the
    :math:`rho` matrix of the graph and its excitations as a propagator
    on the largest eigenvector of the graph. From this, an improved graph
    is stochastically selected, and the process is repeated, lowering
    the energy. If **HDIAG** is specified, it is the hamiltonian matrix
    elements which determine the coupling between determinants, and it
    is the hamiltonian matrix which is diagonalised in each iteration
    in order to calculate the energy.

    .. note:: 
       **GRAPHMORPH** has not been tested with complex wavefunctions.  It will
       almost certainly not work for them.

**VERTEX** **MCDETS**
    Perform Monte Carlo calculations over pure determinant space, which
    is sampled using a series of 'particles' (or 'walkers').

    **MCDETS** is similar to **FCIMC** but maintains at most one
    'particle' at each determinant, which may then contain subparticles
    (which correspond to the individual 'walkers' in **FCIMC**), in
    a binary tree.  This makes some efficiency savings where the same
    information about a determinant is not duplicated.

    **FCIMC** and **MCDETS** calculations share many of the same options
    (see Walker Monte Carlo options, below).

**VERTEX** **RETURNPATHMC**
    Use a spawning algorithm which is constrained in three ways: 

    #. a particle can only be spawned where it will increase its
       excitation level with respect to the reference determinant or
       back to where it was spawned from.
    #. they will spawn back to where their parents were spawned from
       with probability PRet, which is specified using **RETURNBIAS**.
    #. length of spawning chain must be less than the maximum length
       given by **MAXCHAINLENGTH**.

    .. note::
        How can a particle be restricted to spawning to spawning at most
        back to where it was spawned from *and* have a probability of
        spawning back to where its parent was spawed from?
        Documentation *must* be clearer.

    This attempts to circumvent any sign problem in the double
    excitations and the HF, and hopefully this will result in a more stable
    MC algorithm. It remains to be seen if this approach is useful.  Should
    revert to the star graph in the limit of the return bias tending to 1 or
    the length of the spawn chain tending to 1.

    .. note:: 
       **FCIMC**, **GRAPHMORPH**, **MCDETS** and **RETURNPATHMC** have not
       been tested with complex wavefunctions.  It will almost certainly
       not work for them.

       All four are experimental options under development.

Walker Monte Carlo options
--------------------------

The following options are applicable for both the **FCIMC** and **MCDETS** methods:

.. note::
   I have made some guesses on the following option names.  Clearly some keys are broken
   on George's keyboard.  Specifically::

      StepsSft --> STEPSSHIFT
      SftDamp  --> SHIFTDAMP
      DiagSft  --> DIAGSHIFT

   I also had to guess about **BINCANCEL**.  It seems to be a **FCIMC**
   option, but was placed with **MCSTAR** (and was with all the **VERTEX STAR**
   methods).

   This section needs to be extended substantially.

**DIAGSHIFT** [DiagSft]
   Set the initial value of the diagonal shift.

**INITWALKERS** [nWalkers]
    Default 3000.

   Set the initial population of walkers.  
   For CCMC Amplitude, this is the number of samples of the amplitude distribution taken each MC step

**NMCYC** [NMCYC]
   Set the total number of timesteps to take.

**SHIFTDAMP**  [SftDamp]
   Damping factor of the change in shift when it is updated.  <1 means more damping.

**STEPSSHIFT** [StepsSft]
   Default 100.

   Set the number of steps taken before the diagonal shift is updated.

**TAU** [TIMESTEP] 
   Default 0.0.

   For FCIMC, this can be considered the timestep of the simulation. It is a constant which 
   will increase/decrease the rate of spawning/death for a given iteration.

The following options are only available in **FCIMC** calculations:

**READPOPS**
    Read the initial walker configuration from the file POPSFILE.
    **DIAGSHIFT** and **INITWALKERS** given in the input will be
    overwritten with the values read in form POPSFILE.

**SCALEWALKERS** [fScaleWalkers]
    Scale the number of walkers by fScaleWalkers, after having read in data from POPSFILE.

**STARTMP1**
    Set the initial configuration of walkers to be proportional to the MP1 wavefunction. The shift will also
    now be set to the MP2 correlation energy.  This also works in CCMC Amplitude

**GROWMAXFACTOR** [GrowMaxFactor]
    Default 9000.

    Set the factor by which the initial number of particles are allowed to grow before
    they are culled.

**CULLFACTOR** [CullFactor]
    Default 5.

    Set the factor by which the total number of particles is reduced once it reaches the GrowMaxFactor limit

**EQUILSTEPS** [NEquilSteps]
    Default 0
    This indicates the number of cycles which have to
    pass before the energy of the system from the doubles
    population is counted

**SHIFTEQUILSTEPS** [NShiftEquilSteps]
    Default 1000
    This gives the number of iterations after the shift is allowed to change before the shift 
    contributes to the average value printed in column 10.
    The default is 1000 to simply leave out the first few values where the shift is dropping (usually from 0).

**RHOAPP** [RhoApp]
    This is for resummed FCIMC, it indicates the number of propagation steps
    around each subgraph before particles are assigned to the nodes

**SIGNSHIFT**
    This is for FCIMC and involves calculating the change in shift depending on
    the absolute value of the sum of the signs of the walkers.  This should
    hopefully mean that annihilation is implicitly taken into account. Results
    were not too good.

    .. note:: details?  Why "not good"?

**HFRETBIAS** [PRet]
    This is a simple guiding function for FCIMC - if we are at a double
    excitation, then we return to the HF determinant with a probability PRet.
    This is unbiased by the acceptance probability of returning to HF.

    This is not available in the parallel version.

**EXCLUDERANDGUIDE**
    This is an alternative method to unbias for the HFRetBias. It invloves
    disallowing random excitations back to the guiding function (HF
    Determinant).

    This is not available in the parallel version.

**PROJECTE-MP2**
    This will find the energy by projection of the configuration of walkers
    onto the MP1 wavefunction.  DEVELOPMENTAL and possibly not bug-free.

    This is not available in the parallel version.

**FIXPARTICLESIGN**
    This uses a modified hamiltonian, whereby all the positive off-diagonal
    hamiltonian matrix elements are zero. Instead, their diagonals are modified
    to change the on-site death rate. Particles now have a fixed (positive)
    sign which cannot be changed and so no annihilation occurs.  Results were
    not good - this was intended for real-space MC, where large regions of connected
    space were all of the same sign. This is not the case here.
  
    This is not available in the parallel version.

**STARTSINGLEPART**
    This will start the simulation with a single positive particle at the HF,
    and fix the shift at its initial value, until the number of particles gets
    to the INITPARTICLES value.

**MEMORYFACPART** [MemoryFacPart]
    Default 10.D0

    MemoryFacPart is the factor by which space will be made available for extra
    walkers compared to InitWalkers.

**MEMORYFACANNIHIL** [MemoryFacAnnihil]
    Default 10.D0

    MemoryFacAnnihil is a parallel FCIMC option - it is the factor by which space will be 
    made available for annihilation arrays compared to InitWalkers. This generally will need to be 
    larger than memoryfacpart, because the parallel annihilation may not be exactly load-balanced because of 
    non-uniformity in the wavevector and the hashing algorithm. This will tend to want to be larger 
    when it is running on more processors.

**MEMORYFACSPAWN** [MemoryFacSpawn]
    Default 0.5

    A parallel FCIMC option for use with **ROTOANNIHILATION**. This is the factor by which space will be made 
    available for spawned particles each iteration. Several of these arrays are needed for the annihilation 
    process. With **ROTOANNIHILATION**, **MEMORYFACANNIHIL** is redundant, but **MEMORYFACPART** still need to be specified.

**ANNIHILATEONPROCS**
    Default false

    A Parallel FCIMC option. With this, particles are annihilated separately on each node.
    This should mean less annihilation occurs, but it is effectivly running nProcessor
    separate simulations. If there are enough particles, then this should be sufficient.
    Less memory is required, since no hashes need to be stored. Also, no communication is
    needed, so the routine should scale better as the number of walkers grows.

**ROTOANNIHILATION**
    Default false

    A parallel FCIMC option which is a different - and hopefully better scaling - algorithm. 
    This is substantially different to previously. It should involve much less memory.
    **MEMORYFACANNIHIL** is no longer needed (**MEMORYFACPART** still is), and you will need 
    to specify a **MEMORYFACSPAWN** since newly spawned walkers are held on a different array each iteration.
    Since the newly-spawned particles are annihilated initially among themselves, you can still 
    specify **ANNIHILATEATRANGE** as a keyword, which will change things.

**FIXSHELLSHIFT** [ShellFix] [FixShift]
    Default 0,0.D0

    An FCIMC option. With this, the shift is fixed at a value given here, 
    but only for the excitation levels at a value of ShellFix or lower. This will 
    almost definitly give the wrong answers for both the energy and the shift, 
    but may be of use in equilibration steps to maintain particle density at 
    low excitations, before writing out the data and letting the shift change.

**FIXKIISHIFT** FixedKiiCutoff FixShift

    Another fixed shift based approximation method for FCIMC in parallel. However, rather
    than fixing the shift based on an excitation level, it is now fixed according to the 
    Kii value. Determinants lower in energy than FixedKiiCutoff will have their shifts
    fixed to the value given.

**FIXCASSHIFT** [OccCASorbs] [VirtCASorbs] [FixShift]
    Default 0 0 0.D0

    A third fixed shift approximation method for FCIMC in parallel.  In this option, an active
    space is chosen containing a number of highest occupied spin orbitals (OccCASorbs) and a 
    number of lowest unoccupied spin orbitals (VirtCASorbs).  The shift is then fixed (at FixShift)
    for determinants with excitations within this space only.  I.e. determinants for which the spin 
    orbitals lower in energy than the active space are completely occupied and those higher in 
    energy are completely unoccupied.

**SINGLESBIAS** [SinglesBias]
    Default 1.D0

    This represents the factor to which singles are biased towards over double excitations from a determinant.
    This works with the NONUNIFORMRANDEXCITS excitation generators for FCIMC code. Normally, the
    pDoubles is given by the number of doubles divided by the total excitations from HF. Now, 
    the number of singles in the total excitations term is multiplied by SinglesBias. Alternatively,
    SinglesBias can be set to less than 1 to bias towards doubles.

**FINDGROUNDDET**
    Default=.false.

    A parallel FCIMC option. If this is on, then if a determinant is found with an energy lower 
    than the energy of the current reference determinant, the energies are rezeroed and the
    reference changed to the new determinant. For a HF basis, this cannot happen, but with 
    rotated orbital may be important.

**DEFINEDET** [DefDet(NEl)]
    Default=.false.

    A parallel FCIMC option.  This allows the reference determinant to be chosen based on that specified in
    the input with this keyword - rather than the default HF determinant chosen according to the energies of 
    the orbitals.  The determinant is specified by a series of NEl integers (separated by spaces) 
    which represent the occupied spin orbitals.

**DIRECTANNIHILATION**
    Default=.false.

    A parallel FCIMC option. This annihilation algorithm has elements in common with rotoannihilation
    and the default annihilation, but should be faster and better scaling than both of these, with
    respect to the number of processors. There are no explicit loops over processors, and newly-spawned
    particles are sent directly to their respective processors.

**ANNIHILATEEVERY** [iAnnInterval]
    Default=1

    A parallel FCIMC option which can only be used with default annihilation algorihtm. This will
    mean that annihilation will only occur every iAnnInterval iterations.

**ANNIHILATDISTANCE** [Lambda]
    Default=0.D0

    A serial FCIMC option. Here, walkers i and j have the chance to annihilate each other
    as long as they are on connected determinants. They will annihilate with probability
    given by -Lambda*Hij*(Si*Sj). This is hoped to increase annihilation and allow fewer
    particles to be needed to sample the space correctly. When Lambda=0.D0, it should be 
    equivalent to the original annihilation algorithm. Warning - this is much slower than
    normal annihilation.

**ANNIHILATERANGE** [**OFF**]
    Default=.true.

    A parallel FCIMC option. This is a slightly different annihilation algorithm, where only
    one sort of the full set of particles is needed. This should greatly reduce the time needed
    for annihilation of large numbers of particles. However, the load-balancing across processors
    may not be so good. This option is now on by default and can only be switched off via the input
    file by specifying **OFF** after the keyword.

**LOCALANNIHIL** [Lambda]
    
    A parallel FCIMC option. An additional diagonal death rate is included at the annihilation
    stage for particles which are only singly occupied. The probability of death is given by
    Tau*EXP(-Lambda*ExcitDensity) where ExcitDensity is the approximate density of particles in
    the excitation level of the particle. This should raise death through this local annihilation,
    and hence keep the shift at a more resonable value in the undersampled regime. This will
    hopefully mean that a more accurate energy value can be obtained by removing the random
    killing of particles which arises from such a low shift value.

    This is now commented out in the code

**UNBIASPGENINPROJE**
    Default false
    
    An FCIMC serial option. Here, the acceptance probabilities are not unbiased for
    the probability of generating the excitation. Instead, the unbiasing occurs when the 
    walker contributes to the energy estimator. This should reduce the variance for the 
    energy estimator.

**GLOBALSHIFT** **OFF**
    Default true

    This option can only be turned off by specifying **OFF**

    A parallel FCIMC option. It is generally recommended to have this option on. This will 
    calculate the growth rate of the system as a simple ratio of the total walkers on all processors
    before and after update cycle, rather than a weighted average. This however is incompatable with culling, and so 
    is removed for update cycles with this in. This should be more stable than the
    default version and give a more reliable shift estimator for large systems.

**MAGNETIZE** [NoMagDets] [BField]
    Default false

    This is a parallel FCIMC option. It chooses the largest weighted MP1 components and records their 
    sign. If then a particle occupies this determinant and is of the opposite sign, it energy,
    i.e. diagonal matrix element is raised by an energy given by BField. First parameter is an
    integer indicating the number of determinants to 'magnetize', and the second is a real
    giving the amount the energy of a particle should be raised if it is of an opposite sign.
    
**MAGNETIZESYM** [NoMagDets] [BField]
    Default false

    A parallel FCIMC option. Similar to the MAGNETIZE option (same arguments), but in addition to 
    the energy being raised for particles of the opposite sign, the energy is lowered by the same 
    amount for particles of 'parallel' sign.
    
**GRAPHSIZE** [NDets]
    In ResumFCIMC, this is the number of connected determinants to form the
    graph which you take as your sumsystem for the resummed spawning.  Must
    have an associated RhoApp.

**HAPP** [HApp]
    Default 1.

    In ResumFCIMC, this indicates the number of local applications of the
    hamiltonian to random determinants before the walkers are assigned
    according to the resultant vector.

**NOBIRTH**
    Force the off-diagonal :math:`H` matrix elements to become zero,
    and hence obtain an exponential decay of the initial populations
    on the determinants, at a rate which can be exactly calculated and
    compared against. 
    
    This is no longer functional, but commented out in the
    code.

**MCDIFFUSE** [Lambda]
    Default 0.0.

    Set the amount of diffusion compared to spawning in the **FCIMC**
    algorithm.
  
    This is no longer functional and commented out in the code.

**FLIPTAU** [FlipTauCyc]
    Default: off.

    Cause time to be reversed every FlipTauCyc cycles in the **FCIMC**
    algorithm. This might help with undersampling problems.

    This is no longer functional and commented out in the code.

**NON-PARTCONSDIFF**
    Use a seperate partitioning of the diffusion matrices, in which
    the antidiffusion matrix (+ve connections) create a net increase of
    two particles.

    This is no longer functional and commented out in the code.

**FULLUNBIASDIFF**
    Fully unbias for the diffusion process by summing over all connections.

    This is no longer functional and commented out in the code.

**NODALCUTOFF** [NodalCuttoff]
    Constrain a determinant to be of the same sign as the MP1
    wavefunction at that determinant, if the normalised component of
    the MP1 wavefunction is greater than the NodalCutoff value.

    This is no longer functional and commented out in the code.

**NOANNIHIL**
    Remove the annihilation of particles on the same
    determinant step.

**REGENDIAGHELS**
    Default .false.
    This is a parallel FCIMC option, which means that the diagonal hamiltonian matrix
    element for each particle is calculated on the fly, rather than stored with the
    particle. This will free up more memory, but will probably lead to slightly slower
    calculations.

**REGENEXCITGENS**
    This option will regenerate the excitation generator for each particle, every time a 
    new random excitation is wanted. This is MUCH slower for the same number of particles
    (10x?). However, this frees up a lot more memory to store more particles.

**PRINTDOMINANTDETS** [NoDeterminants] [MinExcLevel] [MaxExcLevel]
    Default=.false.

    This is a parallel FCIMC option.  With this keyword, at the end of a calculation a DOMINANTDETS file
    is printed containing the NoDeterminants most populated determinants between excitation
    levels of MinExcLevel and MaxExcLevel (inclusive).  This must be used with rotoannihilation.

**PRINTDOMSPINCOUPLED** [OFF]
    Default=.true.

    This a parallel FCIMC option to go with the one above.  It takes the list of dominant determinants
    chosen based on their populations and adds to the list all the spin coupled determinants that 
    are not already there.  This prevents any spin contamination when we truncate the available 
    determinants.  This is automatically on, but can be turned off using this keyword followed by OFF.

**SPAWNDOMINANTONLY**
    Default=.false.

    This is a parallel FCIMC option.  It takes a DOMINANTDETS file (printed using the above keywords)
    and reads it in at the beginning of the calculation.  During the calculation, if a walker is
    to be spawned with an excitation level of those printed in DOMINANTDETS, this is only allowed if
    the determinant is in the list of dominant determinants.  This does not allow truncation of 
    the doubles, and must be used with rotoannihilation.
    
**STARMINORDETERMINANTS**    
    Default=.false.
    
    This is a parallel FCIMC option.  It goes along with the **SPAWNDOMINANTONLY** keyword.  If this
    is present, spawning to determinants not in the dominant list is done with a star approximation.
    That is, spawning onto minor determinants is allowed, but these walkers are only allowed
    to spawn back to the parent from which they came.  The walkers undergo death and annihilation
    like usual (however, the walkers for annihilation are chosen randomly as they differ depending
    on their parent).

**FINDGUIDINGFUNCTION** [iGuideDets]
    Default=.false. [100]

    This is a parallel FCIMC option.  At the end of a spawning calculation, the iGuideDets most populated
    determinants are found and these and their final walker populations (with sign) are printed out 
    (in order of their bit strings) to a file named GUIDINGFUNC to be used in the subsequent calculation.

**USEGUIDINGFUNCTION** [iInitGuideParts]    
    Default=.false.

    This is a parallel FCIMC option.  This option takes the GUIDINGFUNC file produced in a previous calculation
    and reads in the guiding (or annihilating) function from it.  The population on the HF determinant in this 
    guiding function is then set to be iInitGuideParts, and the remaining determinants are populated based on 
    their occupations from the previous calculation (in GUIDINGFUNC) relative to the HF determinant.
    The function of this guiding function is then to sit in the background of a calculation, able to annihilate
    walkers, but unable to itself spawn of have its walkers die.
    Assuming the GUIDINGFUNC from the previous calculation has the correct nodal structure, this guiding function
    should serve to instantly remove walkers spawned with the incorrect sign.

**TRUNCATECAS** [OccCASOrbs] [VirtCASOrbs]
    This is a parallel FCIMC option, whereby the space will be truncated according to the specified CAS.
    The arguments indicate the active electrons, and then the number of active virtual orbitals.
    These values can be dynamically updated throughout the simulation via use of the CHANGEVARS facility.

**TRUNCINITIATOR** 
    This is a parallel FCIMC option.  The keyword requires an initiator space to first be defined (usually 
    via **TRUNCATECAS**, but could be by **EXCITE**).  
    This is then a variation on a kind of CAS-star approach.  Spawning is subject to the contraint
    that walkers spawned from determinants outside the active space only live if they are being spawned onto 
    determinants that are already occupied.  If walkers spawned on a new determinant have non-initiator parents,
    these spawns are 'aborted'.  A special case is if in the same iteration walkers are spawned on a new 
    determinant both from inside and outside the active space - in this case we treat the active space to have 
    spawned a second earlier, the determinant is then treated as occupied and the non-active space walkers are 
    allowed to live (providing they are the same sign of course).
    NOTE: This is currently only possible using **DIRECTANNIHILATION**.

**DELAYTRUNCINITIATOR** [IterTruncInit]
    This goes with the above.  This allows us to first start with an active space only calculation and then at some
    iteration (given by IterTruncInit), to expand to the **TRUNCINITIATOR** method.  The beginning of the 
    **TRUNCINITIATOR** method may also be started dynamically by putting TRUNCINITIATOR in the CHANGEVARS file.
    At the moment, when this happens, tau is also reduced by a factor of 10.  This should maybe be played with at 
    some stage though.

**KEEPDOUBSPAWNS**
    This keyword goes along with the above **TRUNCINITIATOR**.  This is an extra exception which means that if 
    two determinant spawn on the same determinant with the same sign, they are allowed to live no matter where they
    came from.  This is different from the original case where if both of these had come from non-initiator 
    determinants, they would have both been killed.

**ADDTOINITIATOR** [InitiatorWalkerNo]
    This is again an addition to the above few options.  In this case, if a determinant outside the initiator space
    builds up a significant population (greater than InitiatorWalkerNo), it is treated as being in the initiator
    space and may spawn on occupied or unoccupied determinants as it likes.  This is reassessed at each iteration
    however, so determinant may move in and out of the initiator space as the populations vary.

**INCLDOUBSINITIATOR**
    This keyword also goes with **TRUNCINITIATOR**, and is a parallel FCIMC option.  When it is present, all doubly
    excited determinants are included in the initiator space, and are allowed to spawn as usual.

The following option are only available in **MCSTAR** calculations:

**BINCANCEL** 
    This is a seperate method to cancel down to find the residual
    walkers from a list, involving binning the walkers into their
    determinants. This has to refer to the whole space, and so is
    much slower.  See also the **WAVEVECTORPRINT** and **POPSFILE**
    options in the **LOGGING** block.

**STARORBS** [iStarOrbs] [**NORETURN** | **ALLSPAWNSTARDETS**]
    Default=.false. , NORETURN = OFF

    A parallel FCIMC option. Star orbs means that determinants which 
    contain these orbitals can only be spawned at from the HF determinant, 
    and conversly, can only spawn back at the HF determinant. iStarOrbs is
    the integer variable which decides how many orbitals are in this high-
    energy space, and take the iStarOrbs number of highest energy orbitals
    to construct it. **NORETURN** is an optional keyword specifier. If it
    is specified, then any excitations from the HF to these high-energy
    determinants (doubles) are left to die and cannot respawn back to the
    HF determinant. **ALLSPAWNSTARDETS** is another optional keyword, which
    means that all particles can spawn at determinants with star orbitals, and
    once there, annihilation can occur. However, they cannot respawn anywhere
    else and are left there to die.

**EXCITETRUNCSING** [iHightExcitsSing]
    Default=.false.

    This is a parallel FCIMC option, where excitations between determinants where 
    at least one of the determinants is above iHighExcitsSing will be restricted to be single excitations.

**EXPANDFULLSPACE** [iFullSpaceIter]
    Default=0
    
    This is a parallel FCIMC option. When this is set, the space initially is truncated at excitation level of ICIlevel,
    set by the value of the EXCITE parameter, or the CAS space given by TRUNCATECAS. If EXPANDFULLSPACE is set, then the 
    system will continue to be truncated, but will
    expand to the full space after iteration iFullSpaceIter.
    Hopefully expanding the space in this way will allow quicker
    convergence, without needing to do this dynamically through the use of CHANGEVARS which may be difficult for
    long/queued jobs.

**INITAMPLITUDE** dInitAmplitude

   For CCMC Amplitude, this is the initial amplitude in the Hartree-Fock determinant, and normalization factor for the wavefunction.
   Default 1.0

**CLUSTERSIZEBIAS** dProbSelNewExcitor

   For CCMC Amplitude, this is the probability that the cluster selection algorithm terminates after each addition of an excitor.
   Larger values will bias towards smaller cluster.
   Default 0.7 (Range 0-1)

**NSPAWNINGS** nSpawnings

   For CCMC, this is the number of spawnings attempted from each cluster (unless **EXACTSPAWN** is specified).  Default 1

**SPAWNPROP**
   For Amplitude CCMC use NSPAWNINGS as a total number of spawnings, and distribute them according to the Amplitudes of clusters.
   

Return Path Monte Carlo options
-------------------------------

**MAXCHAINLENGTH** [CLMAX]
    Set the maximum allowed chain length before a particle is forced to
    come back to its origin.

**RETURNBIAS** [PRet]
    Set the bias at any point to spawn at the parent determinant.

Perturbation theory options
---------------------------

**MPTHEORY** [**ONLY**]
    In addition to doing a graph theory calculation, calculate the Moller--Plesset
    energy to the same order as the maximum vertex level from the
    reference determinant (e.g. with 2-vertex sum the MP2 energy is
    obtained, with 3-vertex the MP3 energy etc.  Within the **VERTEX SUM**
    hierarchy, this will only work with **VERTEX SUM HDIAG**.
    In the **VERTEX MC** hierarchy, do a Moller--Plesset calculation 
    instead of a path-integral one.  Requires **HDIAG**, and **BIAS**=0.D0.
    Can be used without a **METHODS** section.  If a **METHODS** section is
    needed to specify different numbers of cycles at each level, then
    **MCDIRECTSUM** must also be set, either in the main block of the **CALC**,
    or by using **VERTEX MCDIRECT** instead of **VERTEX MC**.
    Note that the MP2 energy
    can be obtained in conjunction with a **VERTEX STAR** calculation.

    **ONLY**
        Run only a MP2 calculation.  This is only available when
        compiled in parallel.  The only relevant **CALC** options are the
        **EXCITATIONS** options: all other **CALC** keywords are ignored
        or over-ridden.  No **LOGGING** options are currently applicable.

        Whilst in principle integrals are only used once, this optimal
        algorithm is not currently implemented.  The speed of a **CPMD**-based
        calculation thus benefits from having a **UMatCache** of non-zero size.

        .. warning::
            It is currently assumed that the calculation is restricted.

**EPSTEIN-NESBET**
    Apply Epstein--Nesbet perturbation theory, rather than
    Moller--Plesset.  Only works for **VERTEX SUM NEW** and **VERTEX
    SUM HDIAG** and only at the 2-vertex level.

**LADDER**
    Use ladder diagram perturbation theory, rather than Moller--Plesset.
    The energy denominator is :math:`E_0-E_I+|H_{0I}|^2`.  Only works
    for **VERTEX SUM NEW** and **VERTEX SUM HDIAG** and only at the
    2-vertex level.

**MODMPTHEORY**
    Perform a hybrid of Epstein--Nesbet and Moller--Plesset theory,
    which includes only the :math:`\bra ij||ij ket +\bra ab||ab ket`
    terms in the denominator.  Only works for **VERTEX SUM NEW** and
    **VERTEX SUM HDIAG** and only at the 2-vertex level.

Diagonalisation options
-----------------------

Options for performing a full diagonalisation in the space of the full
basis of spin orbitals.

.. warning::
  This quickly becomes prohibitively expensive as system size increases.

**ACCURACY** [B2L]
    Desired level of accuracy for Lanczos routine.

**BLOCK** [**ON** **OFF**]
    Default off. 

    Determines whether the Hamiltonian is calculated for each block
    or not.  This only works for **COMPLETE**.

**BLOCKS** [NBLK]
    Set number of blocks used in Lanczos diagonalisation.

**COMPLETE**
    Perform a full diagonalisation working out all eigenvectors
    and eigenvalues.  if **HAMILTONIAN** is **OFF**, discard the
    eigenvectors and eigenvalues after having used them for calculation.
    Relevant options are **HAMILTONIAN** and **BLOCK**.

.. note::
  When would it be advantageous to save the eigenvalues and -vectors
  are a diagonalisation?

**EIGENVALUES** [NEVAL]
    Required number of eigenvalues.

**ENERGY**
    Calculate the energy by diagonalising the Hamiltonian matrix.
    Requires one of **COMPLETE**, **LANCZOS**, or **READ** to be set.

    Exact E(Beta) is printed out as:
    
    .. math::
          \text{E(Beta)} = \frac{ \sum_{\alpha} E_{\alpha} e^{-\beta E_{\alpha}} } { \sum_{\alpha} e^{-\beta E_{\alpha}} }

    The result will, of course, change depending upon the symmetry subspace
    chosen for diagonalization for finite temperatures.

    The diagonalization procedure creates a list of determinants, which
    is printed out to the DETS file.

    The weight, :math:`w_{\veci}` and weighted energy, :math:`w_{\veci}
    \tilde{E}_{\veci}` are also calculated for all NPATH determinants.

    .. note::
        **ENERGY** was documented twice in the INPUT_DOC file.  This is not
        particularly helful...  
        
        I have (hopefully) combined them correctly.

**JUSTFINDDETS**

    This is an option to be used in conjunction with **ENERGY** and exact diagonalization methods.
    If specified, the diagonalization routines will just enumerate all the determinants and will
    not try to form the hamiltonian or diagonalize it. No energy will therefore be found, but
    enumerating all the determinants can be useful for histogramming methods in FCIMC methods.

**KRYLOV** [NKRY]
    Set number of Krylov vectors.

**LANCZOS**
    Perform a Lanczos block diagonalisation on the Hamiltonian matrix.  

    Relevant parameters are **BLOCKS**, **KRYLOV**, **ACCURACY**,
    **STEPS** and **EIGENVALUES**.

**READ**
    Read in eigenvectors and eigenvalues of the Hamiltonian matrix from a previous calculation.

**STEPS** [NCYCLE]
   Set the number of steps used in the Lanzcos diagonalisation.

Graph morphing options
----------------------

A new approach developed by George Booth.  Take an initial starting graph
and over many iterations allow the determinants contained within the
graph to change, so that the resultant graph is a better approximation
to the true ground state.

**GRAPHBIAS** [GraphBias]
    If at each iteration the graph is being completely renewed, then this
    bias specifies the probability that an excitation of the previous
    graph is selected to try and be attached, rather than one of the
    determinants in the previous graph.

**GRAPHSIZE** [NDets]
    Specify the number of determinants in the graph to morph.

**GROWGRAPHSEXPO** [GrowGraphsExpo]
    Default is 2.D0. 

    The exponent to which the components of the excitations vector
    and the eigenvector are raised in order to turn them into
    probabilities. The higher the value, the more that larger weighted
    determinants will be favoured, though this might result in the graph
    growing algorithm getting stuck in a region of the space.

**GROWINITGRAPH**
    Grow the initial graph non-stochastically from the excitations of
    consecutive determinants.

**INITSTAR**
    Set up the completely connected two-vertex star graph, and use as
    the starting point for the morphing. 
    
    Automatically changes the NDets parameter to reflect the number of
    double excitations in the system.

**ITERATIONS** [Iters]
    The number of graph morphing iterations to perform.

**MAXEXCIT** [iMaxExcitLevel]
    Limit the size of the excitation space by only allowing excitations
    out to iMaxExcitLevel away from HF reference determinant.

**MCEXCITS** [NoMCExcits]
    Stochastically sample the space of excitations from each determinant in the
    graph with NoMCExcits determinants chosen per determinant.
    For the FCIMC code, this represents the number of attempted spawns per iteration
    in the spawning step.

**MOVEDETS** [NoMoveDets]
    Grow the graphs using an alternative Monte Carlo, where a number
    of determiants are deleted from the previous graph and reattached
    elsewhere in the graph in a stochastic manner, according to the
    probabilities given by the application of the :math:`rho` propagator
    to the eigenvector of the previous graph.

**NOSAMEEXCIT**
    Ignore the connections between determinants which are of the
    same excitation level in comparison to the reference determinant.
    Currently only available in conjunction with **INITSTAR**, so the
    starting graph is simply the doubles star graph (with no cross
    connections).

**ONEEXCITCONN**
    Grow the graph by attaching only determinants which differ by one
    excitation level to the connecting vertex in the previous graph.
    Currently not implemented with MoveDets.

**SINGLESEXCITSPACE**
    Restrict the space into which the current graph is allowed to morph
    to just single excitations of the determinants in the current graph.
    This should reduce the scaling of the algorithm.

Monte Carlo options
-------------------

Options for performing a Monte Carlo calculation on a vertex sum (as
specified in the **METHODS** section).

The Monte Carlo routines have only ever been tested for molucular and
model systems and probably are not currently functional for **CPMD**
or **VASP** based calculations.

See the reports by Ramin Ghorashi ([RGPtIII]_) and George Booth
([GHBCPGS]_).

**CALCVAR**
   Only available for performing full vertex sums using the **HDIAG**
   formulation to evaluate the thermal density matrix elements.

   Calculate a theoretical approximation to the expected variance if a
   non-stochastic MC run were to be performed, with the parameters given,
   at the chosen vertex level.  Currently the expected variance is sent
   to STOUT as a full variance for the total energy ratio.  Causes the
   calculation to take longer since the generation probabilities of
   the graphs must all be calculated.  The sum over graphs of the
   generation probabilities is also printed out for each vertex
   level. This should equal 1, since we are working with normalised
   probabilities.

**POSITION** [IOBS JOBS KOBS]
   Sets the position of the reference particle.

**CIMC**
    Perform a configuration interation space Monte Carlo.

**BETAEQ** [BETAEQ]
    Default is set to be :math:`\beta`, as set above. 

    Set :math:`\beta` to have a different value for the equilibriation steps.

    .. note::
      What are the equilibriation steps?

**BIAS** [G_VMC_FAC]
    Default 16.

    Vertex level bias for **FULL** **MC**. Positive values bias toward
    larger graphs, negative values towards smaller graphs.

    For **SINGLE** and **MULTI** level MC (using a composite 1-vertex
    graph containing a full sum previously performed), this is the
    probability of generating a graph which is not the composite graph.
    The default is invalid, and this must be set manaully.  Stochastic
    time MC is used.  If BIAS is negative, then | BIAS | is used, but
    stochastic-time MC is not performed.

    .. note::
        BIAS seems to do two very different things if it is set to a negative value.
        Please clarify.

**DETSYM** [MDK(I), I=1,4]
    The symmetry of the **CIMC** determinant.

    .. note::
        Specify the symmetry how?

    .. note::    
        If any if the **CIMC** options are set without **CIMC** being
        specified, the code will return an error and exit.**

**EQSTEPS** [IEQSTEPS]
    The number of equilibriation sets for the CI space Monte Carlo routine.

**GRAPHEPSILON** [GRAPHEPSILON]
   Default 0.0.

   The minimum significant value of the weight of a graph.  
   
   Ignore the contributions to the weight and :math:`\tilde{E}` of all
   graphs with a weight that is smaller in magnitude than GRAPHEPSILON.

**IMPORTANCE** [G_VMC_PI] 
    Ddefault 0.95.

    Set the generation probability for the MC routine.  This is the
    probability that new determinants are excitations of the pivot, i.

**MCDIRECTSUM**
   Perform Monte Carlo on graphs summing in energies weighted with the
   weight/generation probability of the graph.


**PGENEPSILON** [PGENEPSILON]
   Default 0.0.

   Set the minimum significant value of the generation probability of a graph. 

   Because for larger graphs, the calculation of the generation
   probability is subject to numerical truncation errors, generation
   probabilities which are lower than a certain value are unreliable,
   and can cause the Monte Carlo algorithm to get stuck: if a graph had a
   very small generation probability, it would be difficult for a Monte
   Carlo run to accept a move to a different graph.  If the magnitude
   of the generation probability of a graph is smaller than PGENEPSILON,
   then a new graph is generated.

   Setting this too high could cause problems in the graph generation phase,
   so NECI will exit with an error if it generates 10000 successive
   graphs each with generation probabilities below PGENEPSILON.

**SEED** [G_VMC_SEED]
    Default -7.

    Set the random seed required for the Monte Carlo calculations.

**STEPS** [IMCSTEPS]
    Set the number of steps for the CI space Monte Carlo routine.

**VVDISALLOW**
   Disallow V-vertex to V'-vertex transitions in stochastic time Monte
   Carlo: i.e. allow only transitions to graphs of the same size.

Weighting schemes
^^^^^^^^^^^^^^^^^

By default the vertex sum Monte Carlo algorithm selects excitations
with no bias.  The variance of a Monte Carlo calculation can be reduced
by preferentially selecting for certin types of excitation.

**EXCITWEIGHTING** [g_VMC_ExcitFromWeight g_VMC_ExcitToWeight G_VMC_EXCITWEIGHT] [g_VMC_ExcitToWeight2]
   Default 0.d0 (unweighted) for all values.

   A weighting factor for the generation of random excitations in the
   vertex sum Monte Carlo.  A parameter set to zero has a corresponding
   weighting factor of 1.

   For generating an excitation from occupied spin orbitals i and j to
   unoccupied spin orbitals k and l:

       * the probability of choosing pair (ij) is proportional to 
           .. math::
                e^{(E_i+E_j) \text{g\_VMC\_ExcitFromWeight} }.

       * the probability of choosing pair (kl) is proportional to 
           .. math::
                e^{-(E_k+E_l) \text{g\_VMC\_ExcitToWeight}} e^{|\bra ij|U|kl\ket|*\text{G\_VMC\_EXCITWEIGHT}} |E_i+E_j-E_k-E_l|^{\text{g\_VMC\_ExcitToWeight2}}.

**POLYEXCITWEIGHT** [g_VMC_ExcitFromWeight g_VMC_PolyExcitToWeight1 g_VMC_PolyExcitToWeight2 G_VMC_EXCITWEIGHT]
    Default 0.0 for all values (i.e. unweighted: all weighting factors
    are set to 1).

    Weighting system for the choice of virtual orbitals in
    the excitations.

    The probability of choosing the pair of spin orbitals, kl, to excite
    to is set to be constant for :math:`E_k+E_l` is less than
    g_VMC_PolyExcitToWeight1.  For higher energy virtual orbitals,
    the weighting applied is a decaying polynomial which goes as:

        .. math::
           (E_k+E_l-\text{g\_VMC\_PolyExcitToWeight1}+1)^{-\text{g\_VMC\_PolyExcitToWeight2}}

    g_VMC_PolyExcitToWeight1 is constrained to be not more than the
    energy of the highest virtual orbital.

**POLYEXCITBOTH** [g_VMC_PolyExcitFromWeight1 g_VMC_PolyExcitFromWeight2 g_VMC_PolyExcitToWeight1 g_VMC_PolyExcitToWeight2 G_VMC_EXCITWEIGHT]
    Identical to **POLYEXCITWEIGHT**, except that the polynomial weighting
    function applies also to the occupied orbitals.  This means that there
    is another variable, since now the 'ExcitFrom' calculation also needs
    a value for sigma, and for the exponent.  The sigma variables are
    now both under similar constraints as specified above, which means
    that they cannot be larger or smaller than the highest and lowest
    energy orbital respectivly.  This prevents the PRECALC block from
    getting stuck, or from finding local variance minima.

    .. note::
         What is sigma?

**CHEMPOTWEIGHTING** [g_VMC_PolyExcitFromWeight2 g_VMC_PolyExcitToWeight2 G_VMC_EXCITWEIGHT]
    Weighting is of the same form as POLYEXCITBOTH, but sigma is
    now constrained to be at the chemical potential of the molecule.
    Has only two parameters with which to minimise the expected variance.

**CHEMPOT-TWOFROM** [g_VMC_ExcitWeights(1) g_VMC_ExcitWeights(2) g_VMC_ExcitWeights(3) G_VMC_EXCITWEIGHT]
    When choosing the electron to excite, use a a increasing polynomial
    up to the chemical potential and a decaying polynomial for spin orbitals
    above the chemical potential, in order to encourage mixing of
    the configurations around the HF state. Contains three adjustable
    parameters and testing needs to be done to see if this is
    beneficial. Expected to make more of a difference as the vertex
    level increases.

   .. note::
         What is the actual weighting form of **CHEMPOT-TWOFROM**?

**UFORM-POWER**
    New power form for the U-matrix element weighting using the
    appropriate **EXCITWEIGHT** element, which is believed to be
    better. This uses the form :math:`W=1+|U|^{\text{EXCITWEIGHT}}`, rather than the
    exponential form.

**STEPEXCITWEIGHTING** [g_VMC_ExcitWeights(1) g_VMC_ExcitWeights(2) G_VMC_EXCITWEIGHT]
    This excitation weighting consists of a step function between the HF virtual and occupied electon manifold (i.e. step is at the chemical potential)
    When choosing an electron to move, the weight for selecting the electron is increased by 1 if the electron oribital has energy above the chemical potential
    and by g_VMC_ExcitWeights(1,1) if below. This occurs for both electrons. When choosing where to excite to, the situation is reversed, and the weight of selecting the
    unoccupied orbital is increased by 1 if the orbital is a hole in the occupied manifold and g_VMC_ExcitWeights(2,1) if a virtual orbital in the occupied manifold. 
    Bear in mind that the parameters are NOT probabilities. If we are at a higher excitation level w.r.t. HF, then more electrons will be in the virtual manifold, 
    which will alter the normalisation, and mean that when selecting electrons to excite, there will be an increasingly small probability of selecting them from the 
    occupied manifold. The opposite is true when choosing where to put them.

    Simply put, if the parameters are both < 1, then the biasing will preferentially generate excitations which reduce the excitation level.
    
    U-weighting is the third parameter as before.



Experimental options
--------------------

.. note::
  More documentation on these options needed.

**EXCITATIONS** **FORCEROOT**
   Force all excitations in **VERTEX** [**SUM** **STAR**] **NEW**
   calculations to come from the root.

**EXCITATIONS** **FORCETREE**   
   Disallow any excitations in a **VERTEX** **SUM** **NEW** which are
   connected to another in the graph, forcing a tree to be produced.
   Not all trees are produced however.

**FULLDIAGTRIPS**
    An option when creating a star of triples, to do a
    full diagonalisation of the triples stars, without any
    prediagonalisation. Very very slow...


**LINEPOINTSSTAR** [LinePoints]
    Set the number of excited stars whose eigenvalues are evaluated when 
    using StarStars, in order to determine linear scaling.

**NOTRIPLES**
    Disallow triple-excitations of the root determinant as the 3rd vertex
    in **HDIAG** calculations at the third vertex level and higher.
