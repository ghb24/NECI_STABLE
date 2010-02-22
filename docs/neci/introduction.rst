.. _introduction:

============
Introduction
============

NECI is a rapidly developing code based on a post Hartree--Fock electronic structure method.

It calculates electron correlation via path-resummations in Slater Determinant
space [SumPaper]_, [StarPaper]_, [ThomPhDThesis]_.  More recent work has
focused on a novel Quantum Monte Carlo theory, FCIQMC, that obtains the Full
Configuration Interaction energy.  See [FCIQMC]_ and [Initiator]_ for more details.

As a standalone package, NECI can perform calculations on electrons confined to
a box, the uniform electron gas and the hubbard model.  

NECI can also read in wavefunctions, or a set of integrals based on them, of
molecular systems produced by another program (e.g. [DALTON]_ or [MolPro]_) and
run calculations using them as the basis for forming the necessary Slater
Determinants.

Finally, NECI can also be compiled as a library for integration into existing
codes.  Currently this has been performed for the [CPMD]_ or [VASP]_ plane-wave
packages, allowing calculations to be performed on periodic systems.
