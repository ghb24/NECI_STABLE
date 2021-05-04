title: LOGGING Block 
---

### LOGGING Block

The LOGGING block specifies the output of the calculation and which
status information of the calculation shall be collected. This block is
started with the `logging` keyword and terminated with the `endlog`
keyword.

-   <span style="color: blue">**logging**  
    </span> Starts the LOGGING block.

-   <span style="color: blue">**endlog**  
    </span> Terminates the LOGGING block.

-   <span style="color: blue">**hdf5-popsfile**  
    </span> Sets the format to read and write the wave function to HDF5.
    Requires building with the `ENABLE-HDF5` cmake option.

-   **popsfile \(n\)**  
    Save the current wave function on disk at the end of the
    calculation. Can be used to initialize subsequent calculations and
    continue the run. This is enabled by default. \(n\) is optional and,
    when given, specifies that every \(n\) iteration, the wave function
    shall be saved. Setting \(n=-1\) disables this option.

-   **popsFileTimer \(n\)**  
    Write out a the wave function to disk every \(n\) minutes, each time
    overwriting the last output.

-   **hdf5-pops-write**  
    Sets the format to write the wave function to HDF5. Requires
    building with the `ENABLE-HDF5` cmake option.

-   **hdf5-pops-read**  
    Sets the format to read the wave function to HDF5. Requires building
    with the `ENABLE-HDF5` cmake option.

-   **highlyPopWrite \(n\)**  
    Print out the \(n\) most populated determinants at the end of the
    calculation. Is enabled by default with \(n=15\).

-   **replicas-popwrite [\(n\)]**  
    Have each replica print out its own highly populated determinants in
    a separate block instead of one collective output of the on average
    most important determinants. Optionally, \(n\) is the numbers of
    determinants to print, this is the same \(n\) as in
    `highlypopwrite`.

-   **inits-exlvl-write \(n\)**  
    Sets the excitation level up to which the number of initiators is
    logged to \(n\). Defaults to \(n=8\).

-   **binarypops**  
    Sets the format to write the wave function to binary.

-   **nomcoutput**  
    Suppress the printing of iteration information to stdout. This data
    is still written to disk.

-   **stepsOutput \(n\)**  
    Write the iteration data to stdout/disk every \(n\) iterations.
    Defaults to the number of iterations per shift cycle. Setting
    \(n=0\) disables iteration data output to stdout and uses the shift
    cycle for disk output.

-   **fval-energy-hist**  
    Create a histogram of the scaling factor used for the auto-adaptive
    shift over the energy of a determinant. Only has an effect if
    `auto-adaptive-shift` is used.

-   **fval-pops-hist**  
    Create a histogram of the scaling factor used for the auto-adaptive
    shift over the population of a determinant. Only has an effect if
    `auto-adaptive-shift` is used.

#### Semi-stochastic output options

-   **write-core**  
    When performing a semi-stochastic calculation, adding this option to
    the Logging block will cause the core space determinants to be
    written to a file called CORESPACE. These can then be further read
    in and used in subsequent semi-stochastic options using the
    `read-core` option in the CALC block.

-   **write-most-pop-core-end \(n\)**  
    At the end of a calculation, output the \(n\) most populated
    determinants to a file called CORESPACE. This can further be read in
    and used as the core space in subsequent calculations using the
    `read-core` option.

#### RDM output options

These options control how the RDMs are printed. For a description of how
the RDMs are calculated and the content of the files, please see section
<a href="#sec:rdms" data-reference-type="ref" data-reference="sec:rdms">1.7</a>.

-   **calcRdmOnfly \(i\) \(step\) \(start\)**  
    Calculate RDMs stochastically over the course of the calculation.
    Starts sampling RDMs after \(start\) iterations, and outputs an
    average every \(step\) iterations. \(i\) indicates whether only
    1-RDMs (1), only 2-RDMs (2) or both are produced.

-   **rdmLinSpace \(start\) \(n\) \(step\)**  
    A more user friendly version of `calcrdmonfly` and
    `rdmsamplingiters`, this samples both 1- and 2-RDMs starting at
    iteration \(start\), outputting an average every \(step\) iterations
    \(n\) times, then ending the calculation.

-   **diagFlyOneRdm**  
    Diagonalise the 1-RDMs, yielding the occupation numbers of the
    natural orbitals.

-   **printOneRdm**  
    Always output the 1-RDMs to a file, regardless of which RDMs are
    calculated. May compute the 1-RDMs from the 2-RDMs.

-   **writeRdmsToRead off**  
    The presence of this keyword overrides the default. If the OFF word
    is present, the unnormalised `TwoRDM_POPS_a***` files will
    definitely not be printed, otherwise they definitely will be,
    regardless of the state of the `popsfile/binarypops` keywords.

-   **readRdms**  
    This keyword tells the calculation to read in the `TwoRDM_POPS_a***`
    files from a previous calculation. The restarted calc then continues
    to fill these RDMs from the very first iteration regardless of the
    value put with the `calcRdmOnFly` keyword. The calculation will
    crash if one of the `TwoRDM_POPS_a***` files are missing. If the
    `readRdms` keyword is present, but the calc is doing a
    `StartSinglePart` run, the `TwoRDM_POPS_a***` files will be ignored.

-   **noNormRdms**  
    This will prevent the final, normalised `TwoRDM_a***` matrices from
    being printed. These files can be quite large, so if the calculation
    is definitely not going to be converged, this keyword may be useful.

-   **writeRdmsEvery \(iter\)**  
    This will write the normalised `TwoRDM_a***` matrices every \(iter\)
    iterations while the RDMs are being filled. At the moment, this must
    be a multiple of the frequency with which the energy is calculated.
    The files will be labelled with incrementing values -
    `TwoRDM_a***.1` is the first, and then next `TwoRDM_a***.2` etc.

-   **write-spin-free-rdm**  
    Output the spin-free 2-RDMs to disk at the end of the calculation.

-   **printRoDump**  
    Output the integrals of the natural orbitals to a file.

-   **print-molcas-rdms**  
    It is now possible to calculate stochastic spin-free RDMs with the
    GUGA implementation. This keyword is necessary if one intends to use
    this feature in conjunction with `Molcas` to perform a spin-free
    Stochastic-CASSCF. It produces the three files `DMAT, PAMAT` and
    `PSMAT`, which are read-in by `Molcas`.
