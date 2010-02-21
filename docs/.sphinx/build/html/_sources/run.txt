.. _run:

---
Run
---

NECI
----

.. note::
    How to obtain FCIDUMP/density-fitting input files?

.. code-block:: bash

    $ neci.x input_file

If no file is given, then it takes input options from STDIN.  This is rarely useful, however.

NECI prints output to STDOUT, so output needs to be captured in some way:

.. code-block:: bash

    $ neci.x input_file > output_file
    $ neci.x nput_file | tee output_file

CPMD-NECI
---------

The converged Kohn--Sham orbitals obtained from a **OPTIMIZE
WAVEFUNCTION** CPMD calculation can be used as input for a NECI
calculation.

In contrast to the molecular case, NECI calculations based upon
CPMD-generated wavefunctions are called from within CPMD itself.
This allows us to take advantage of many routines that CPMD already
possesses (FFT routines, initialisation, reading in the wavefunctions
etc.).

To run, specify **QMC** in the **&CPMD** section of the CPMD input file.
**RESTART WAVEFUNCTIONS OCCUPATION DENSITY COORDINATES LATEST** must
also be specified.  Running CPMD (assuming it has been correctly compiled
with the appropriate NECI library) then calls NECI to read the NECI
input file and perform the desired calculation.  

For gamma-point calculations:

.. code-block:: bash

    $ gcpmd.x input_file > output_file

For k-point calculations:

.. code-block:: bash

    $ kcpmd.x input_file > output_file

There are many other appropriate options that can be specified in the
CPMD input file rather than the NECI input file.  Please see the CPMD
manula and the local CPMD documentation detailing addtions the Alavi
group has made.

Interacting with running calculations
-------------------------------------

It is useful to be able to interact with a running calculation, particularly
during iterative and cyclic processes and is not available (nor suitable) for
all types of calculation.  Currently this functionality is implemented for
**FCIMC** calculations.

NECI checks for a file called CHANGEVARS in the working directory of the
calculation.  CHANGEVARS can be placed on any node when run on multi-node
compute servers.  The CHANGEVARS file is read and echoed to STDOUT and
then deleted so that it does not affect any subsequent calculations.

Valid options to CHANGEVARS are:

**EXCITE** [excitation_level]
    Change the maximum excitation level of determinants included in the
    simulation.  A value less than or equal to 0 or greater than the number of
    electrons in the system sets the maximum excitation level to use the full
    determinant space.
**TRUNCATECAS** [OccCASOrbs] [VirtCASOrbs]
    Change the space used to the specified complete active space (CAS). A value
    equal to or less than 0 or greater than the number of electrons in the system
    sets it to the full determinant space.
**SOFTEXIT**
    Stop calculation as soon as possible.  This facility is useful when running
    on machines with fixed walltimes, especially when used with the watchdog.py
    script in the utils/ subdirectory.
**WRITEPOPS**
    Print out the current walker population to POPSFILE.
**VARYSHIFT**
    Exit fixed shift phase and allow the shift to vary according to [FCIQMC]_.
**NMCYC** [ncycles]
    Change the number of Monte Carlo cycles to perform.
**TAU** [tau]
    Change the timestep, tau, for the simulation.
**DIAGSHIFT** [shift]
    Change the shift.
**SHIFTDAMP** [damping]
    Change the damping parameter used to adjust the shift.
**STEPSSHIFT** [nsteps]
    Change the number of Monte Carlo cycles performed between updating the shift.
**SINGLESBIAS** [bias]
    Change the bias for generating single excitations over double excitations
    when using the non-uniform random excitation generator.
**ZEROPROJE**
    Rezero the averaged energy estimators.  This is useful when the initial value of
    the energy estimators are a long way from the instantaneous values, causing a 
    slow convergence of the averaged values.
**ZEROHIST**
    Rezero the averaged histogramming vectors.
**PARTIALLYFREEZE** [nPartFrozen nHolesFrozen]
    Change the maximum number of holes, nHolesFrozen, allowed in the nPartFrozen
    number of spin-orbitals in the core valence region.
    Determinants with a larger number of holes in the lowest nPartFrozen
    spin-orbitals are not considered.
    See the input option in the **INTEGRALS** section for more details.
**PARTIALLYFREEZEVIRT** [nVirtPartFrozen nElVirtFrozen]
    Similar to **PARTIALLYFREEZE**, allow only Slater determinants with
    at most nElVirtFrozen electrons in the nVirtPartFrozen number of virtual
    spin-orbitals.
**PRINTERRORBLOCKING**
    Print the blocking analysis.
**STARTERRORBLOCKING**
    Start the blocking analysis.
**RESTARTERRORBLOCKING**
    Restart the blocking analysis.
**PRINTSHIFTBLOCKING**
    Print the shift blocking analysis.
**RESTARTSHIFTBLOCKING**
    Restart the shift blocking analysis.
**EQUILSTEPS** [ncycles]
    Change the number of initial Monte Carlo cycles to ignore in the averaging
    of the energy and the shift.
**STARTHIST**
    Begin histogramming the determinant populations if the tCalcFCIMCPsi
    is on and the histogramming has been set up.
**HISTEQUILSTEPS** [ncycles]
    Change the iteration at which the histogramming begins to the value
    specified.
**TRUNCINITIATOR**
    Expand the CAS calculation to a **TRUNCINITIATOR** calculation if
    **DELAYTRUNCINITIATOR** is present in the input.
**ADDTOINIT** [nwalkers]
    Will change the cutt-off population for which walkers are added to the
    initiator space.  The population must be above specified value.

Many of these options are also valid options in the main input file and are covered in
more depth in input.
