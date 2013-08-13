.. _input_logging:

-------
Logging
-------

.. note::
 All the logging options should say to which file they print output.  Please correct this!

**LOGGING**
  Start the logging input block.  The logging options allow additional
  (potentially expensive, potentially verbose) information to be
  printed out during a calculation.  By default, all logging options
  are turned off.

[Logging options---see below.]

**ENDLOG**
    End the logging input block.

General options
---------------

**NOMCOUTPUT**
    
    Ensures no output to stdout from the fcimc or ccmc iterations
    
**FMCPR** [**LABEL**, **RHO**, **1000**, **EXCITATION**]
    More than one of the options can be specified.

    Log the following to the PATHS file:

    **LABEL**
       Logs the determinants contained by each graph as each determinant
       is generated in the format:
       [:math:`(D_0),(D_1),...,D_v),`]
       where each determinant given as a comma-separated list of the
       indices of the occupied orbitals:
       e.g. :math:`D_0 =` (    1,    2,    9,   10,).

       If CSFs are being used, then the CSF is printed.  There is no newline after this.

       For **MULTI MC** or **SINGLE MC**, only the non-composite graphs are printed.

    **EXCITATION**
       Log each graph in excitation format instead of full format above.  The format is
            [A(    i,    j)->(    a,    b),B(    k,    l)->(    c,    d),...,C(    m,    n)->(    e,    f)]
       where 
            A, B,..., C are determinants in the graph from which the excitation is made. 
            i, j,... are the orbitals within that determinant which are excited from, where (i<j, k<l,...).
            a, b,... are the orbitals they are excited to, where (a<b, c<d, ...).

       This format does not in general provide a unique way of
       specifying multiply connected graphs, but the first possible
       determinant to which the next det in the graph is connected is
       chosen, so what is output should be unique.  Single excitations
       are written as e.g. A(    i,    0)->(    a,    0).

    **RHO**
       Log the :math:`\rho` matrix for each graph in the form: 
           (:math:`\rho_{11}, \rho_{12}, \cdots, \rho_{1v},| \rho_{21}, \rho_{22}, \cdots, \rho_{2v},| \cdots | \rho_{v1}, \rho_{v2}, \cdots, \rho_{vv},|`), 

       where the graph consists of :math:`v`  vertices.  A newline is appended.

    **XIJ**
       Log the xij matrix, which contains the generation probabilities
       of one determinant in the graph from all the others.  For MC this
       is already generated, but for full sums this must be generated,
       so will be slower.  Generation probabilities are set with the
       **EXCITWEIGHTING** option.  

       The format is:
           {:math:`x_{11}, x_{12}, ..., x_{1v},| ... | x_{v1}, x_{v2}, ..., x_{vv},|`}

       In general :math:`x_{ij} \ne x_{ji}`.  The :math:`x_{kk}` element lists
       the number of possible excitations from :math:`k` determinant.
       The matrix is followed by a newline.

       After all these possible options, the following are printed:
            Weight [pGen] ETilde*Weight Class [Accepted]

       pGen is only printed for: Monte Carlo, or if doing a full sum
       and the XIJ logging option is set.  Accepted is only printed
       for Monte Carlo calculations.  A newline is placed at the end
       of this data.  For Monte Calo calculations, the values printed
       depend on the options.  If **LABEL** is set, then all generated
       graphs and their values are printed, otherwise only values of
       accepted graphs are printed.

**CALCPATH** [**LABEL** **RHO**]
    Log CALCPATH_R to PATHS.  Either just label logging or also
    log the :math:`\rho` matrix, with the same format as above.

**HAMILTONIAN**      
    Log HAMIL, PSI and PSI_LONG.

**HFBASIS**
    Log HFBASIS.

**HFLOGLEVEL** [LEVEL]
   Default 0.  
   
   If LEVEL is set to be positive, the density matrices, fock matrices and
   eigenvectors during a Hartree--Fock calculation are printed out to SDOUT.

**MCPATHS**     
    Log MCPATHS data to the MCPATHS file for full vertex sum and MCSUMMARY
    file when using a METHODS section.  Also log to the RHOPII file.

**PSI**
    Log PSI_COMP.

**TIMING** [iGlobalTimerLevel | **LEVEL** iGlobalTimerLevel | **PRINT** nPrintTimer]
   **LEVEL** iGlobalTimerLevel
       Default 40.
       Timing information is only recorded for routines with level less than
       or equal to iGlobalTimerLevel.  Less than 10 means general high level
       subroutines. Greater than 50 means very low level.  Routines without
       a level are always timed (which is most of them).  The greater the value
       of iGlobalTimerLevel, the more routines are timed.  This can affect 
       performance in some cases.
   **PRINT** nPrintTimer
       Default 10.
       Print out timing information for the nPrintTimer routines which took the longest time.

**XIJ**
   Synonym for **FMCPR XIJ**.

**DETS**
**DETERMINANTS**
   Log the list of determinants to DETS and SymDETS if they're calculated.

FCIMC options
-------------

**HISTSPAWN** [iWriteHistEvery]

    This option will histogram the spawned wavevector, averaged over all previous iterations. 
    It scales horrifically and can only be done for small systems which can be diagonalized. 
    It requires a enumeration of all determinants initially to work. It can write out the 
    average wavevector every iWriteHistEvery.
    If a diagonalization option is set, SymDets will also be written out, containing the exact 
    wavevector in the same format from the diagonalization.

**HISTPARTENERGIES** [BinRange] [iNoBins] [OffDiagBinRange] [OffDiagMax]

    This is a histogramming option. It is slow, so not for use unless the diagnostic is needed. It will histogram
    the diagonal hamiltonian matrix element for three types of particle. Two input values are needed. The first
    argument is a real value to give the width of the histogram bin. The second is the number of bins needed (integer).
    Three histograms are produced: EVERYENERGYHIST - this is the histogram over all iterations of every particle in the
    system. ATTEMPTENERGYHIST - this is the histogram of the energy of all attempted spawned particles (including the 
    ones which are successfully spawned). For this one, the contibution to the energy is actually 1/Prob of generating. 
    SPAWNENERGYHIST - this is the histogram of all successfully spawned particles. All these histograms are normalized to
    one before printing out.
    Also now, the off-diagonal matrix elements are histogrammed. OffDiagBinRange is a real input parameter which indicates
    the range of the bins and OffDiagMax is the maximum matrix element to histogram. The doubles and singles will be done
    seperately, as are the accepted spawns and total spawns. Therefore, four files are produced - SINGLESHIST, ATTEMPTSINGLESHIST,
    DOUBLESHIST, ATTEMPTDOUBLESHIST. Again, these are normalized and the ATTEMPT files histogram proportionally to 1/probability
    of generating the excitation.

**AUTOCORR** [NoACDets(2)] [NoACDets(3)] [NoACDets(4)]
    This is a parallel FCIMC option. It will output the histogrammed occupation number for certain
    determinants every iteration. This is so that a seperate standalone ACF program can be used on it.
    Currently the histogramming is evaluated for the HF determinants by default, but can also 
    histogram determinants from other excitation levels. Firstly, it will calculate the 'NoACDets(2)' 
    largest-weighted MP1 components (double excitations). It will then take the largest weighted double
    and do a new MP1 calculation with it as the root. It will then histogram the 'NoACDets(3)' largest 
    weighted triple excitations, and the 'NoACDets(4)' largest quadruple excitations from this calculation
    to also histogram.

**REDUCEDPOPSFILE** [iWritePopsEvery] [iPopsPartEvery]
    This works in the same way as the normal popsfile, but only every iPopsPartEvery particle is printed out.

**POPSFILE** [iWritePopsEvery]
    Default: on.  iWritePopsEvery is an optional argument.

    Write out the necessary information to restart a calculation, including the
    population of the walkers on each determinant.  This is written to 
    POPSFILE or (if **BINARYPOPS** is specified) POPSFILEHEAD and POPSFILEBIN.

    If iWritePopsEvery is supplied, then the determinant populations are
    printed out every iWritePopsEvery Monte-Carlo cycles.  iWritePopsEvery
    should idealy be a multiple of **STEPSSHIFT**, the number of cycles between
    updates to the diagonal shift performed in the **FCIMC** calculation, to
    make sure the start of the next simulation follows smoothly.

    A calculation can then be restarted at a later date by reading the
    determinants back in using **READPOPS** in the **CALC** section. 
    Walker number can also be scaled up/down by using **SCALEWALKERS**.
    If the iWritePopsEvery argument is negative, then the POPSFILE is never
    written out, even at the end of a simulation. This is useful for very large
    calculations where the POPSFILE will take a long time to write out and use
    a lot of disk space.

**BINARYPOPS**
    This means that the popsfile (full or reduced) will now be written out in binary format. 
    This should now take up less disk space, and be written quicker. It can be read in as
    normal without specifying any extra criteria. Two files will be produced, a formatted
    file with the header info and a POPSFILEBIN with the walker information.

**INCREMENTPOPS**
    Append a unique suffix to the POPSFILE* restart file(s) to avoid
    overwriting them.  Note that this can quickly fill up hard drives if used
    with **POPSFILE** iWritePopsEvery: use with care!

**POPSFILETIMER** [PopsfileTime]
    Write out a POPSFILE every 'PopsfileTime' hours of the calculation. Can be used with
    **INCREMENTPOPS** to save previous files.

**ZEROPROJE**
    This is for FCIMC when reading in from a POPSFILE. If this is on, then the energy 
    estimator will be restarted.

**WAVEVECTORPRINT**
    This is for Star FCIMC only - if on, it will calculate the exact eigenvector and
    values initially, and then print out the running wavevector every
    WavevectorPrint MC steps. However, this is slower.

**PRINTFCIMCPSI**
    This works for parallel FCIMC. This will enumerate all excitations (up to the truncation level specified,
    or the full space if not specified), and then histogram the spawning run, writing out the final
    averaged wavefunction at the end.

**HISTEQUILSTEPS** [NHistEquilSteps]
    Default=.false. [0]
    This works when the evolving wavefunction is to be histogrammed (for example using the above 
    **PRINTFCIMCPSI** option, or the **USECINATORBS** orbital rotation option).  
    This sets the histogramming to only begin after NHistEquilSteps iterations.  This is so that the 
    fluctuation populations at the beginning of a calculation may be left out.

**PRINTORBOCCS**
    Default=.false.
    This turns on the histogramming of the determinant populations, and at the end of the spawning, calls a 
    routine to add up the contribution of each orbital to the final wavefunction.  A ORBOCCUPATIONS file is then
    printed containing the orbitals and their normalised absolute occupations.

**WRITEDETE** [NoHistBins] [MaxHistE]
    This is an FCIMC option and will write out a histogram of the energies of determinants which have
    had particles spawned at them and their excitation level. The histogram logs the total
    amount of time spent at a determinant and its energy for each energy range. This is diagnostic 
    information. The first variable to input is the number of histogram bins which will be calculated,
    and the second is the maximum determinant energy of the histogram.

**PRINTTRICONNECTIONS** [TriConMax] [NoTriConBins]
    This is a parallel FCIMC option. It looks at sets of connected determinants i,j and k.  A sign coherent
    triangular connection is one where walkers spawned all around the triangle return to the original 
    determinant with the same sign.  Sign incoherent connections are those where the sign is reversed.  
    If this option is on, two files are printed.  TriConnTotal monitors the number of sign coherent and sign 
    incoherent triangles over the course of the simulation, as well as the sum of the Hij x Hik x Hjk values,
    and the ratios for each.  (The ratios are coherent / incoherent).  TriConnHist prints out a histogram of
    the Hij x Hik x Hjk values for coherent (col 1 and 2) and incoherent (col 3 and 4) triangles.  The histogram
    goes from 0 -> +/- TriConMax with NoTriConBins for each.

**HISTTRICONNELEMENTS** [TriConHElSingMax] [TriConHElDoubMax] [NoTriConHElBins]
    This option histograms all the H elements involved in the triangular connections of determinants mentioned
    above.  These are separated into doubles and singles, and an extra file, containing only the Hjk elements is
    also included.
    The histogram range is between +/-TriConHElSingMax for the singles and +/-TriConHElDoubMax for the doubles, with
    NoTriConHElBins bins for each.
    With this option, some stats are also printed in the output regarding the average magnitudes for each type of H
    elements.

**PRINTHELACCEPTSTATS**
    This option prints out a file (HElsAcceptance) containing information about the nature of the H elements 
    resulting in accepted and not accepted spawns.  This includes the number of not accepted spawns vs accepted, 
    and the average size of the H element involved in accepted and not accepted spawns. 

**PRINTSPINCOUPHELS**
    Default=.false.
    When attempting to spawn on a determinant i, this option finds the determinant j which is spin coupled to i, and 
    prints out a set of stats relating to the sign and magnitude of the H element connecting i and j, Hij.
    These stats are printed in a file named SpinCoupHEl.

**BLOCKEVERYITER**
    Default=.false.
    This will block the projected energy every iteration with the aim of achieving accurate error estimates. 
    Two caveats - it does not take into account the serial correlation between the numerator and denominator of the energy
    expression, and does require a small amount of additional communication each iteration.

**CCMCDEBUG** iCCMCDebug
    Specify the CCMC debug level.  Default 0 (no debugging information printed).  Higher numbers will generate more
    information.

**FCIMCDEBUG** iFCIMCDebug
    Specify the FCIMC debug level.  Default 0 (no debugging information printed).  Higher numbers will generate more
    information.

**CCMCLOGTRANSITIONS** [**NONUNIQUE** **UNIQUE**]
    Do we log all transitions in CCMC.  Very slow and memory intensive - only possible for extremely small systems.
      Default is **UNIQUE**.  If **NONUNIQUE** is specified, then clusters with different orders are distinguished.

GraphMorph options
------------------

**DISTRIBS**
    Write out the distribution of the excitations in each graph as it
    morphs over the iterations. The first column is the iteration number, and
    then subsequent columns denote the number of n-fold excitations in
    the graph.

PRECALC options
---------------

**PREVAR**
    Print the vertex level, Iteration number, parameter, and expected
    variance, for each parameter which was searched for in the **PRECALC**
    block, showing the convergence on the optimum value, to the PRECALC
    file.

**SAVEPRECALCLOGGING**
   Allows different logging levels to be used in the **PRECALC** block
   than for the main calculation.

   All logging options specified before **SAVEPRECALCLOGGING** are only
   used in the the **PRECALC** part of the calculation.  All logging
   options specified after  **SAVEPRECALCLOGGING** are only used in the
   the main part of the calculation.

Monte Carlo options
-------------------

**BLOCKING**
    Perform a blocking analysis on the MC run.  An MCBLOCKS file will be
    produced, which lists log(2)[blocksize], the average of the blocks,
    the error in the blocks(where the blocks are the energy ratio),
    and the full error, treating the energy estimator as a correlated
    ratio of two quantities.

**ERRORBLOCKING** [OFF]
    Default= ErrorBlocking.true. 
    This can be used to turn off the error blocking analysis that is peformed 
    by default on parallel FCIMC calculations.  The default error blocking 
    begins when the sum of the HF population over an update cycle reaches 1000.  
    At the end of the simulation a BLOCKINGANALYSIS file is printed containing
    a list of block sizes with the resulting average of the projected energies 
    calculated over an update cycle, the error in this energy and the error on 
    the calculated error due to the block size.

**BLOCKINGSTARTHFPOP** [HFPopStartBlocking]
    Default=1000
    This can be used to change the HF population that triggers the start of the
    error blocking analysis.  Using this keyword over rides the default, and 
    the blocking starts when the sum of the HF pop over an update cycle reaches 
    HFPopStartBlocking.

**BLOCKINGSTARTITER** [IterStartBlocking]
    Default=.false.
    This can be used to set the error blocking to begin at iteration number 
    IterStartBlocking, rather than a particular HF population.

    The error blocking may also be initiated instantly by using **STARTERRORBLOCKING**
    in the CHANGEVARS file.  Additionally, **PRINTERRORBLOCKING** will print the 
    BLOCKINGANALYSIS file at that point, yet the calculation (and blocking) will 
    continue (note - this file will be overwritten when the calculation ends and the 
    final blocking stats are printed, so it must be renamed if it is to be kept).  
    **RESTARTERRORBLOCKING** in the CHANGEVARS file zeroes all the 
    blocking arrays and starts again from that point in the calculation.

**SHIFTERRORBLOCKING** [OFF]    
    Default= ShiftErrorblocking.true.
    This can be used to turn off the default error blocking of the shift values.  
    This only starts when the shift begins to vary, and may be restarted or the 
    current SHIFTBLOCKINGANALYSIS file printed at that point using CHANGEVARS.

**SHIFTBLOCKINGSTARTITER** [IterShiftBlock]
    This can be used to specify the number of iterations after the shift is allowed
    to change that the shift error blocking begins.

**VERTEX** [**EVERY** n]
    Log the vertex MC with :math:`\tilde{E}` every n (real) cycles
    and/or log the vertex MC contribution every cycle.  Setting
    Delta :math:`=\tilde{E}-\tilde{E}_{\textrm{ref}}`, where
    :math:`\tilde{E}_{\textrm{ref}}` is usually the 1-vertex graph:

    **EVERY**
        write a VMC file with the following info, with a new line each
        time the current graph changes:

             tot # virt steps, # steps in this graph, #verts, Class, Weight, Delta, <sign(W)>, <Delta sign(W)>, ~standard deviation <Delta sign>/<sign>,pgen 
    n:
        write a VERTEXMC file with the following info:

            0, #graphs, <sign(W)>, stdev(sign(W)), <Delta>, <sign Delta>/<sign>, <Delta^2>, acc ratio, trees ratio, nontree+ ratio, non-tree- ratio, <Delta sign(W)>, E~ reference, #sequences,w reference

.. note::
 George, what are most of these values?

**WAVEVECTORPRINT** [nWavevectorPrint]
    Relevant only for Monte Carlo star calculations.
    
    Calculate the exact eigen-vectors and -values initially, and 
    print out the running wavevector every nWavevectorPrint Monte Carlo
    steps. This is slows the calculation down substantially.


Rotate Orbs Options
-------------------

**ROFCIDUMP** [OFF]
    At the end of an orbital rotation (or in the case of a softexit), by default 
    a ROFCIDUMP file will be printed using the transformation coefficients.
    This may then be read in to a spawning calculation.
    In the case of ROFCIDUMP OFF, no FCIDUMP will be printed.
    Note: When reading in the ROFCIDUMP, the number of electrons must be reduced 
    by the number frozen in the previous rotation, and the number frozen set to 0.

**ROHISTOGRAMALL**
    If this keyword is present, two files are printed for all possible histograms.
    One labelled HistHF*, and one HistRot* containing the histogram before and after rotation.
    With this, certain histograms may be turned off by using the below keywords.  
    Alternatively combinations of the keywords below may be used to just print a selection
    of the possible histograms.

**ROHISTOFFDIAG** [OFF]
    Histograms <ij|kl> terms before and after rotation where i<k and j<l.

**ROHISTDOUBEXC** [OFF]
    Histograms the 2<ij|kl>-<ij|lk> terms, the off diagonal hamiltonian elements for double 
    excitations.

**ROHISTSINGEXC** [OFF]
    Histograms the single excitation hamiltonian elements.

**ROHISTER** [OFF]
    Histograms the <ii|ii> values before and after rotation.

**ROHISTONEElINTS** [OFF]
    Histograms the one electron integral terms, <i|h|i>.

**ROHISTONEPARTORBEN** [OFF]
    Histograms the one particle orbital energies, epsilon_i = <i|h|i> + sum_j [<ij||ij>],
    where j is over the occupied orbitals only.

**ROHISTVIRTCOULOMB** [OFF]
    Histograms the two electron coulomb integrals <ij|ij> where i and j are both virtual spatial orbitals
    and i<j.
    
**TRUNCROFCIDUMP** [NoFrozenOrbs]    
    This option goes along with the **USEMP2VDM** rotation option.  Having diagonalised the MP2VDM
    matrix to get the transformation matrix.  This option truncates the virtual orbital space by removing
    the NoFrozenOrbs SPIN orbitals with the lowest occupation numbers (MP2VDM eigenvalues).  Only the 
    remaining orbitals are transformed and included in the ROFCIDUMP that is printed.
    This kind of transformation requires different ordering of the orbitals to that which is standard for 
    spawning calculation, so it is not possible to go straight from this rotation into a spawning calc.
    The ROFCIDUMP must be printed out then read back in.

**WRITETRANSFORMMAT** 
    Default false.
    This keyword must be included if we are doing a natural orbital rotation, and we want to print out
    an MOTRANSFORM file.  This file contains the transformation matrix in binary which can be used with 
    Qchem to get the cube files for the new orbitals.  NOTE: This file is only printed correctly if NECI
    is compiled using PGI when the file is printed. 

Reduced Density Matrix (RDM) Options
-------------------------------------

Currently the 2-RDMs can only be calculated for closed shell systems.  However, calculation and 
diagonalisation of only the 1-RDM is set up for either open shell or closed shell systems.

The theory behind the calculation of the RDMs can be found in my (Deidre's) thesis (Chapter 7).

The only thing that differs significantly in the more recent code is that, for efficiency, the 
diagonal elements of the RDMs (and explicit connections to the HF determinant) 
are only calculated during the iteration whereby the energy or RDM itself is required.  
This is either every time the energy is printed, or at the very end of the calculation.  
Therefore, for an expensive calculation, it is worth only 
calculating the energy at the end, or only a few times throughout the calculation, to avoid the N^2 operation 
required for each determinant to fill the diagonal RDM elements.

The calculation of the diagonal elements is done by keeping track of the average walker populations of each 
occupied determinant, and how long it has been occupied.  The diagonal element from Di is then calculated 
as <Ni> x <Ni> x [No. of iterations Di has been occupied], and this is included every time a determinant becomes 
unoccupied, at the end of the calculation, or if the energy from the RDMs is to be calculated. 

Because the average occupation is accumulated while the RDMs are being calculated (and is not currently zero-ed 
unless a determinant becomes unoccupied), calculations with the energy calculated at different frequencies will 
have very slightly different RDMs.  
This difference appears to be too small to be a problem, but is noted here to avoid unnecessary debugging. 

The average populations are also used for the stochastic elements, which are accumulated throughout the 
calculation of the RDMs.  Further explanation of how and why this is done can be found in my thesis.

**CALCRDMONFLY** [RDMExcitLevel] [RDMIterStart] [RDMEnergyIter]
    This is the main keyword for calculating the RDMs from an FCIQMC wavefunction.  It requires 3 integers.
    The first refers to the type of RDMs to calculate.  A value of 1 will calculate only the 1-RDM.  Any other 
    value will calculate the 2-RDM (which contains the information of the 1-RDM).  The second integer 
    is the number of iterations after the shift has begun to change that we want to begin filling the RDMs.  
    Finally, if the 2-RDMs are being calculated, the RDM energy will be automatically obtained at the 
    end of the calculation.  The 3rd integer refers to how often (every RDMEnergyIter iterations) we want 
    to additionally calculate and print the energy during the calculation.  This will be ignored if only 
    calculating the 1-RDM. 
    Clearly making RDMEnergyIter very large will mean the energy is only calculated with a softexit, or this can 
    also be achieved by using **CALCRDMENERGY** OFF.

    The RDM energy is one measure of the accuracy of the RDMs.  Also printed by default are the maximum error in the 
    hermiticity (2-RDM(i,j;a,b) - 2-RDM(a,b;i,j)) and the sum of the absolute errors.  

Types of RDM calculations

    Using the above keyword, a stochastic RDM calculation on the entire space will be performed, but with 
    the single and double connections to the HF included explicitly.  The type of calculation can be 
    changed by including any of the following keywords.

**EXPLICITALLRDM**
    This performs a completely explicit calculation of the RDMs.  It considers all single and double excitations 
    of each determinant and therefore adds in every occupied connection at each iteration.  It is very 
    expensive and can only be done for very small systems.  Cannot use **HPHF** with this type of calculation.

    If **HISTSPAWN** is also present in the Logging block, the wavefunction will be histogrammed from the same 
    iteration we begin to fill the RDMs, and the RDMs will be constructed using these histogrammed coefficients.

**HFREFRDMEXPLICIT**
    This effectively calculates the matrix leading to the projected energy.  It considers the HF as a reference 
    and explicitly considers all connections to it (in one direction only - so the matrix is not hermitian).  This 
    will clearly not give a variational energy.  If this RDM was calculated using the instantaneous occupations 
    (rather than the occupied average), the energy printed at every iteration would be the same as Eproj.

**HFSDRDM**
    This calculates the RDM amongst the HF and single and double excitations only.  The connections to the HF will 
    be considered explicitly, but connections between singles and doubles stochastically.  This is hermitian, and 
    should give a variational energy.  Cannot use **HPHF** with this type of calculation.

**HFSDREFRDM**
    This effectively calculates a multireference version of the projected energy, using the HF, singles and doubles 
    as a reference.  Again the connections to the HF are explicit, but all others (up to 4-fold excitation) are 
    included stochastically.  Like **HFREFRDMEXPLICIT**, this matrix will not be hermitian and the energy not variational.  
    Cannot use **HPHF** with this type of calculation.

**RDMGHOSTCHILD** REMOVED 
    This option is discussed in my thesis, but has been removed because it is virtually never used, doesn't actually 
    work very well for most systems, and complicated the code a bit.  The version where it is removed was noted in the 
    logs if it needs to be put back in.

Options referring to the 1-RDM.

**DIAGFLYONERDM**
    This option can be used when calculating either the 1- or 2-RDM.  If we're calculating the 2-RDM, the 1-RDM is 
    constructed and then diagonalised (to get the natural orbital occupation numbers - NO_OCC_NUMBERS, and 
    transformation matrix - NO_TRANSFORM).  This cannot be used with either of the non-hermitian options 
    **HFREFRDMEXPLICIT** or **HFSDREFRDM**.
    If this keyword is present the correlation entropy is also calculated and printed in the output.  

**NONOTRANSFORM**
    This option is used if we want to diagonalise the 1-RDM to get the correlation entropy, but don't want to 
    print the NO_TRANSFORM matrix.

**PRINTRODUMP**
    If this keyword is present, the natural orbital transformation matrix will be used to transform the 4-index 
    integrals etc, to produce a new FCIDUMP (ROFCIDUMP) in the natural orbital basis.
    This is quite slow and expensive and probably wants to be avoided unless actually necessary.
    Options for truncating the orbitals based on NO occupation number has been removed from this code, as well 
    as rotating virtual and occupied orbitals separately.  On the todo list is to put this stuff back in and 
    merge these routines with the old equivalents (in RotateOrbs and NatOrb).
    Note: If you print this from a frozen core calculation, the ROFCIDUMP will be printed as though the
    frozen electrons don't exist.  To restart in the rotated basis, you need to unfreeze the core, and reduce 
    the number of electrons in the input by the number originally frozen.

**PRINTONERDM**
    This means the 1-RDM will be constructed and printed, even if we are only really calculating the 2-RDM.

**DUMPFORCESINFO**
    This constructs the symmetry-packed arrays of density matrices required for molpro force calculations, along
    with the Lagrangian term, printed in binary to fciqmc_forces_info.  The readable version is also currently written
    to the end of the OUTPUT file.  We also print out information on the 'hermiticity error' in the Lagrangian
    which we explicitely symmetrise.

Reading in / Writing out the RDMs for restarting calculations.
    
    Two types of 2-RDMs can be printed out.  The final normalised hermitian 2-RDMs of the form TwoRDM_a***, or the 
    binary files TwoRDM_POPS_a***, which are the unnormalised RDMs, before hermiticity has been enforced.  The 
    first are the matrices to be used for F12 calculations etc (these 2-RDM(i,j;a,b) matrices are printed in 
    spatial orbitals with i<j, a<b and i,j<a,b).  The second are the ones to read back in if a calculation 
    is restarted (they are also printed in spatial orbitals with i<j and a<b, but for both i,j,a,b and a,b,i,j 
    because they are not yet hermitian).  These are the matrices exactly as they are at that point in the calculation.  
    By default the final normalised 2-RDMs will always be printed, and the TwoRDM_POPS_a*** files are connected to the 
    POPSFILE/BINARYPOPS keywords - i.e. if a wavefunction POPSFILE is being printed and the RDMs are being filled, 
    a RDM POPSFILE will be also.
    If only the 1-RDM is being calculated, OneRDM_POPS/OneRDM files will be printed in the same way.
    The following options can override/modify these defaults.

**WRITERDMSTOREAD** [OFF]
    The presence of this keyword overrides the default.  If the OFF word is present, the unnormalised TwoRDM_POPS_a*** 
    files will definitely not be printed, otherwise they definitely will be, regardless of the state of the 
    POPSFILE/BINARYPOPS keywords.  

**READRDMS**
    This keyword tells the calculation to read in the TwoRDM_POPS_a*** files from a previous calculation.  The 
    restarted calc then continues to fill these RDMs from the very first iteration regardless of the value put with 
    the **CALCRDMONFLY** keyword.  The calculation will crash if one of the TwoRDM_POPS_a*** files are missing.  If 
    thisn **READRDMS** keyword is present, but the calc is doing a **STARTSINGLEPART** run, the TwoRDM_POPS_a*** files 
    will be ignored.

**NONORMRDMS**
    This will prevent the final, normalised TwoRDM_a*** matrices from being printed.  These files can be quite 
    large, so if the calculation is definitely not going to be converged, this keyword may be useful.

**WRITERDMSEVERY** [IterWriteRDMs]
    This will write the normalised TwoRDM_a*** matrices every IterWriteRDMs iterations while the RDMs are being 
    filled.  At the moment, this must be a multiple of the frequency with which the energy is calculated.  The 
    files will be labelled with incrementing values - TwoRDM_a***.1 is the first, and then next TwoRDM_a***.2 etc.
    This option currently only works for the 2-RDMs.

