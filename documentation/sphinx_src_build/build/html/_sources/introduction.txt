.. _introduction:

============
Introduction
============

NECI is a rapidly developing code based on a post Hartree--Fock electronic structure method.

It calculates electron correlation via path-resummations in Slater Determinant space [SumPaper]_, [StarPaper]_, [ThomPhDThesis]_.

As a standalone package, NECI can perform calculations on electrons confined to a box, the uniform electron gas and the hubbard model.  

NECI can also read in wavefunctions, or a set of integrals based on them, of molecular systems produced by another program (e.g. [DALTON]_ or [MolPro]_) and run calculations using them as the basis for forming the necessary Slater Determinants.

Finally, NECI can also be compiled as a library for integration into existing codes.  Currently this has been performed for the [CPMD]_ or [VASP]_ plane-wave packages, allowing calculations to be performed on periodic systems.



.. [SumPaper] A combinatorial approach to the electron correlation problem, Alex J. W. Thom and Ali Alavi, J. Chem. Phys. 123, 204106, (2005).
.. [StarPaper] Electron correlation from path resummations: the double-excitation star, Alex J. W. Thom, George H. Booth, and Ali Alavi, Phys. Chem. Chem. Phys., 10, 652-657 (2008).
.. [ThomPhDThesis]  Towards a quantum Monte Carlo approach based on path resummations, Alex J.W. Thom, PhD Thesis (2006).
.. [DALTON] DALTON, a molecular electronic structure program, Release 2.0 (2005), see http://daltonprogram.org/
.. [MolPro] MOLPRO is a package of ab initio programs written by H.-J. Werner, P. J. Knowles, R. Lindh, F. R. Manby,  M. Schütz, P. Celani, T. Korona, G. Rauhut, R. D. Amos, A. Bernhardsson, A. Berning, D. L. Cooper, M. J. O. Deegan, A. J. Dobbyn, F. Eckert, C. Hampel, G. Hetzer, A. W. Lloyd, S. J. McNicholas, W. Meyer, M. E. Mura, A. Nicklaß, P. Palmieri, R. Pitzer, U. Schumann, H. Stoll, A. J. Stone, R. Tarroni, and T. Thorsteinsson, see http://www.molpro.net
.. [CPMD] CPMD, http://www.cpmd.org/, Copyright IBM Corp 1990-2008, Copyright MPI für Festkörperforschung Stuttgart 1997-2001.
.. [VASP] VASP
