title: INTEGRAL Block 
---
### INTEGRAL Block

The INTEGRAL block can be used to freeze orbitals and set properties of
the integrals. The block is started with the `integral` keyword and
terminated with the `endint` keyword.

-   **integral**  
    Starts the INTEGRAL block.

-   **endint**  
    Terminates the INTEGRAL block.

-   **freeze \(n\) \(m\)**  
    Freeze \(n\) core and \(m\) virtual orbitals which are not to be
    considered active in this calculation. The orbitals are selected
    according to orbital energy, the \(n\) lowest and \(m\) highest
    orbitals in energy are frozen.

-   **freezeInner \(n\) \(m\)**  
    Freeze \(n\) core and \(m\) virtual orbitals which are not to be
    considered active in this calculation. The orbitals are selected
    according to orbital energy, the \(n\) highest and \(m\) lowest
    orbitals in energy are frozen.

-   **partiallyFreeze \(n_\text{orb}\) \(n_\text{holes}\)**  
    Freeze \(n_\text{orb}\) core orbitals partially. This means at most
    \(n_\text{holes}\) holes are now allowed in these orbitals.

-   **partiallyFreezeVirt \(n_\text{orb}\) \(n_\text{els}\)**  
    Freeze \(n_\text{orb}\) virtual orbitals partially. This means at
    most \(n_\text{els}\) electrons are now allowed in these orbitals.

-   **hdf5-integrals**  
    Read the 3-body integrals for a transcorrelated ab-initio
    Hamiltonian from an HDF5 file. Ignored when the
    `molecular-transcorrelated` keyword is not given. Requires compiling
    with HDF5.

-   **sparse-lmat**  
    Store the 3-body integrals in a sparse format to save memory.
    Initialisation and iterations might be slower. Requires
    `hdf5-integrals`.
