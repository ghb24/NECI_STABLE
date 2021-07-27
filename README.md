# Codebase for NECI.

FCIQMC code developed by George Booth and Ali Alavi, 2013
Copyright (c) 2013, Ali Alavi
Please see "LICENSE" file for GNU GPL v.3 license associated with this software.

This program is integrated in Molpro with the permission of
George Booth and Ali Alavi.

Please feel free to contact myself (George Booth) or Ali Alavi for
queries or help at george.booth@kcl.ac.uk and office-alavi@fkf.mpg.de.

Contributers to the Alavi group NECI codebase, in alphabetical order;

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


