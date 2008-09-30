.. _run:

---
Run
---

NECI
----

.. note::
  How to obtain FCIDUMP/density-fitting input files?

.. code-block:: bash

 neci.x input_file

If no file is given, then it takes input options from STDIN.  This is rarely useful, however.

NECI prints output to SDTOUT, so output needs to be captured in some way:

.. code-block:: bash

 neci.x input_file > output_file
 neci.x nput_file | tee output_file


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

 gcpmd.x input_file > output_file

For k-point calculations:

.. code-block:: bash

 kcpmd.x input_file > output_file

There are many other appropriate options that can be specified in the
CPMD input file rather than the NECI input file.  Please see the CPMD
manula and the local CPMD documentation detailing addtions the Alavi
group has made.

Soft exits
----------

Soft exits allow a calculation to be halted prematurely yet still
perform any necessary post-processing of the calculation and print out,
for example, the memory and timing information.  NECI checks the working
directory of the calculation for a file called SOFTEXIT, which can be
created by the user by, for instance, using touch:

.. code_block:: bash

   touch SOFTEXIT

On creation of SOFTEXIT, the next time SOFTEXIT is checked to see
if it exists, the code will exit cleanly and quietly, with no error
messages.  The SOFTEXIT file is deleted so that it does not affect 
any subsequent calculations.

This functionality is especially suited to iterative and cyclic
processes and is not available (nor suitable) for all types of
calculation.

Currently, SOFTEXIT applies only to the following calculation
type(s):

* **FciMC** 
  If SOFTEXIT is created, then the current iteration is completed,
  all post-processing of the calculation up to the current iteration
  is performed.
