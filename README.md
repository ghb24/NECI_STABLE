# Codebase for NECI.

This is the codebase of NECI,
a state-of-the-art implementation of the Full Configuration Interaction Quantum Monte Carlo (FCIQMC) algorithm,
a method based on a stochastic application of the Hamiltonian matrix on a sparse sampling of the wave function.
The program utilizes a
very powerful parallelization and scales efficiently to more than 24 000 central processing unit cores.

## Contributors
Contributors to the Alavi group NECI codebase, in alphabetical order

Ali Alavi, Robert Anderson, Nick Blunt, George Booth, Deirdre Cleland, Werner Dobrautz,
Khaldoon Ghanem, Kai Guther, Peter Jeszenski, Niklas Liebermann, Florian Merz,
Jennifer Mohr, Catherine Overy, Simon Smart, James Shepherd, James Spencer, Anthony Stone,
Lauretta Schwarz, Alex Thom, Robert Thomas, David Thompson, Oskar Weser.

## Documentation

The documentation can be found [here](https://www2.fkf.mpg.de/alavi/neci/devel/index.html).

## Compilation

To compile the code, clone the repository (now referred to as `source_directory`)
and create a `build_directory` outside the `source_directory`.

Then you can type
```bash
cd ${build_directory}
cmake ${source_directory}
make -j
```
and are ready to go.
`cmake` will try to find all dependencies and fails if it cannot.

For detailed compilation and run options, please see the user documentation.


## Interfaces

Since NECI is mainly a program to solve the CI-problem
and only depends on the 1- and 2-electronic integrals.
Other code is required to compute these integrals and perform
other optimization tasks (like SCF orbital rotations).

So far it has been interfaced to
[OpenMolcas](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rasscf.html#stochastic-casscf-method),
[Molpro](https://www.molpro.net/),
[PySCF](https://pyscf.org/), and [VASP](https://www.vasp.at/).


## License

FCIQMC code developed by George Booth and Ali Alavi, 2013
Copyright (c) 2013, Ali Alavi
Please see "LICENSE" file for GNU GPL v.3 license associated with this software.

## Questions

Please feel free to contact myself (George Booth) or Ali Alavi for
queries or help at george.booth@kcl.ac.uk and office-alavi@fkf.mpg.de.
