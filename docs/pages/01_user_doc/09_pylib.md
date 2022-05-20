---
title: Using the NECI pylib
---

## Using the NECI pylib

A python frontend is available in the `neci_guga` python library, which
is based on the `neci` build target (i.e. plain neci, without complex
or multi-replica support) and is built by executing
```bash
  make neci_guga_pylib
```
in the neci build directory. This will create a python3 library
`neci_guga.<build-specifier>.so`, which is installed in the
`python/` subdirectory of the neci build directory.

To use it, load the `python` subdirectory of the neci build directory
into the library path of python, either by
```bash
  export PYTHONPATH=<neci_build>/python:$PYTHONPATH
```
or by adding
```bash
  import sys
  sys.path.append('<neci_build>/python')
```
to the calling python module or script.

The `neci_guga` python module can then be loaded in calling python code
with `import neci_guga` and provides the following functionality:

-   `neci_guga.init_guga(S)`
    Takes the desired total spin `S` and initializes the GUGA functionality
    of neci by reading in an existing FCIDUMP file.

-   `neci_guga.clear_guga()`
    Clears all memory and deletes all objects initialized by `init_guga`,
    returns 0 on success, 1 else.

-   `neci_guga.csf_matel(D_i, D_j)`
    Returns the matrix element between `D_i` and `D_j`, passed as an array
    of the size of the number of electrons in the DefineDet format.

-   `neci_guga.diag_matel(D)`
    Returns the diagonal matrix element of a CSF `D`, passed as an array
    of the size of the number of electrons in the DefineDet format.

-   `neci_guga.run_neci(perm)`
    Reads a neci input file `neci.inp` and an FCIDUMP file in the current
    directory, using an orbital permutation `perm` to re-order the orbitals
    used for the calculation (with respect to the FCIDUMP file). Then, a neci
    calculation with the specified input is run and the weight of the leading
    CSF is returned.
    The permutation is given by specifying the new position for each orbital,
    i.e. a permutation
    ```python
      perm = [3, 4, 2, 1]
    ```
    would put orbital 1 in the third position, orbital 2 in the fourth position,
    orbital 3 in the second position and orbital 4 in the first position.
