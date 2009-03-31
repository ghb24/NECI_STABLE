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

FCIMC options
-------------

**HISTDETS** [iWriteHistEvery]

    This option will histogram the spawned wavevector, averaged over all previous iterations. 
    It scales horrifically and can only be done for small systems which can be diagonalized. 
    It requires a diagonalization initially to work. It can write out the average wavevector every iWriteHistEvery.
    If done, SymDets will also be written out, containing the exact wavevector in the same format from the 
    diagonalization.

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
    Default: on.  Default iWritePopsEvery (optional argument) 100000.
    Print out the determinants every iWritePopsEvery Monte-Carlo cycles.
    iWritePopsEvery should idealy be a multiple of **STEPSSHIFT**, the number of 
    cycles between updates to the diagonal shift performed in the 
    **FCIMC** calculation, to make sure the start of the next simulation follows
    smoothly

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

**ZEROPROJE**
    This is for FCIMC when reading in from a POPSFILE. If this is on, then the energy 
    estimator will be restarted.

**WAVEVECTORPRINT**
    This is for Star FCIMC only - if on, it will calculate the exact eigenvector and
    values initially, and then print out the running wavevector every
    WavevectorPrint MC steps. However, this is slower.

**WRITEDETE** [NoHistBins] [MaxHistE]
    This is an FCIMC option and will write out a histogram of the energies of determinants which have
    had particles spawned at them and their excitation level. The histogram logs the total
    amount of time spent at a determinant and its energy for each energy range. This is diagnostic 
    information. The first variable to input is the number of histogram bins which will be calculated,
    and the second is the maximum determinant energy of the histogram.



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
    
