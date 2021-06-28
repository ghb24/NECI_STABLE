---
title: Davidson RAS code
---

## Davidson RAS code

NECI has an option to find the ground state of a RAS space using a
direct CI davidson approach, which does not require the Hamiltonian to
be stored. This code is particularly efficient for FCI and CAS spaces,
but is less efficient for CI spaces.

To perform a davidson calculation, put

<!-- NOTE putting bash for NECI commands, at least for now -->
```bash
davidson ras1 ras2 ras3 ras4 ras5
```

in the Methods block, inside the Calc block. The parameters ras1-ras5
define the RAS space that will be used. These are defined as follows.
First, split all of the spatial orbitals into theree sets, RAS1, RAS2
and RAS3, so that RAS1 contains the lowest energy orbitals, and RAS3 the
highest. Then, ras1, ras2 and ras3 define the the number of spatial
orbitals in RAS1, RAS2 and RAS3. ras4 defines the minimum number of
electrons in RAS1. ras4 defines the maximum number of electrons in RAS3.
These 5 parameters define the ras space.

This method will allocate space for up to 25 Krylov vectors. It will
iterate until the norm of the residual vector is less than \(10^{-7}\).
If this is not achieved in 25 iterations, the calculation will simply
stop and output whatever the current best estimate at the ground state
is.

This code should be able to perform FCI or CAS calculations for spaces
up to around \(5\times10^6\) or so, but will probably struggle for
spaces much larger than this.

The method has only been implemented with RHF calculations and with
\(M_s=0\).
