<img src="./docs/media/neci.jpg" align="left" alt="NECI logo" width="150">

# Codebase for NECI.

This is the codebase of NECI,
a state-of-the-art implementation of the Full Configuration Interaction Quantum Monte Carlo (FCIQMC) algorithm,
a method based on a stochastic application of the Hamiltonian matrix on a sparse sampling of the wave function.
The program utilizes a
very powerful parallelization and scales efficiently to more than 24 000 central processing unit cores.

## Contributors
Contributors to the NECI codebase, in alphabetical order

Ali Alavi, Robert Anderson, Nick Blunt, George Booth, Deirdre Cleland, Werner Dobrautz,
Khaldoon Ghanem, Kai Guther, Philip Haupt, Peter Jeszenski, Giovanni Li Manni,
Niklas Liebermann, Florian Merz, Jennifer Mohr, Catherine Overy, Arta Safari, Pradipta Samanta,
Simon Smart, James Shepherd, James Spencer, Anthony Stone, Lauretta Schwarz, Pierre-Louis Taillat,
Alex Thom, Robert Thomas, David Thompson, Eugenio Vitale, Oskar Weser.

## Documentation

The documentation can be found [here](https://www2.fkf.mpg.de/alavi/neci/stable/index.html).

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
and only depends on the 1- and 2-electronic integrals
other code is required to compute these integrals and perform
other optimization tasks (like SCF orbital rotations).

So far it has been interfaced to
[OpenMolcas](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rasscf.html#stochastic-casscf-method),
[Molpro](https://www.molpro.net/),
[PySCF](https://pyscf.org/), and [VASP](https://www.vasp.at/).


## License

FCIQMC code initially developed by George Booth and Ali Alavi, 2013
Copyright (c) 2013, Ali Alavi
Please see "LICENSE" file for GNU GPL v.3 license associated with this software.

## Questions

Please feel free to create an issue at the public GitHub page
if you have any questions.

## Acknowledgments


| | |
| --- | --- |
| <img src="./docs/media/MPG_logo.png" alt="MPG logo" align="left" width="400">  | Currently most of the development of `NECI` happens at the [Max Planck Institute for Solid State Research](https://www.fkf.mpg.de/de).  |
| <img src="./docs/media/cambridge_university2.svg" alt="Cambridge logo" align="left" width="300">  | The early development of `NECI` was performed at the University of Cambridge. |
| <img src="./docs/media/UKRI_logo.png" alt="ESPRC logo" align="left" width="400"> |  The UK [EPSRC](https://www.ukri.org/councils/epsrc/) Grants (EP/I014624/1 and EP/J003867/1) supported the early development of `NECI`. |
| <img src="https://github.com/zulip/zulip/blob/main/static/images/logo/zulip-icon-512x512.png" alt="Zulip logo" align="center" width="50"> | [Zulip](https://zulip.com/) is an open-source modern team chat app designed to keep both live and asynchronous conversations organized. They support this project with free access to Zulip Cloud Standard. |
