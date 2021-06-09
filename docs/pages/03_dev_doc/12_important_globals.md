---
title: Guide to specific code in NECI
---

# Guide to specific code in NECI

## Important modules and common (global) variable names

### Inexplicable names and anachronisms

A number of the variable names in NECI are largely inexplicable, or
confusing, outside of their origin in historical accident. This section
tries to clarify some of these.

-   **ARR, BRR**<br>
    `ARR(:,1)` contains a list of spin orbital (Fock) energies in order
    of increasing energy.

    `ARR(:,1)` contains a list of spin orbital (Fock) energies indexed
    by the spin-orbital indices used in the calculation.

    `BRR` contains a list of spin orbital indices in increasing (Fock)
    energy order. Where multiple degenerate orbitals have the same
    symmetry, they are clustered so that spin orbitals with the same
    symmetry are adjacent to each other, and within that ordered by
    \(m_s\) value.

-   **TotWalkers and TotParts**<br>
    `TotWalkers` refers to the number of determinants or sites that must
    be looped over in the main list. This may include a number of blank
    slots if the hashed storage is being used.

    The number of particles (walkers) on a core is stored in `TotParts`.
    The total number of particles in the system (including all nodes) is
    stored in `AllTotParts`.

-   **InitWalkers**<br>
    The number of particles *per processor* before a simulation will
    enter variable shift mode.

-   **G1**<br>
    The variable name `G1` is a historical anachronism. It is an array
    containing symmetry information about given basis functions. In
    particular `G1(orb)` contains information about the one electron
    orbital `orb`.

    The data is of type `BasisFN`:

    ```Fortran
    type :: symmetry
        sequence
        integer(ints64) :: S
    end type
    type :: BasisFn
        type(symmetry) :: sym
            !! Spatial symmetry
        integer :: k(3)
            !! K-vector
        integer :: Ms
            !! Spin of electron. Represented as +/- 1
        integer :: Ml
            !! Magnetic quantum number (Lz)
        integer :: dummy
            !! Padding for alignment
    end type
    ```

    Not all of these values are required for all simulations. The `sym`
    element points at a padded 64-bit integer. This is done for memory
    alignment reasons that are no longer important.

## Parallelism

FCIQMC is a highly parallelisable algorithm. Implementationally this is
achieved through the use of independent MPI processes that communicate
using MPI.

The raw `MPI_*` routines, provided by MPI should not be used directly
inside NECI for a couple of reasons:

Depending on the build configuration and compiler, many variables may
change their size or alignment, leading to code which only works on some
compilers, and

Naming conventions for linking vary between compilers.

The `Parallel_neci` module abstracts these implementational details
away, presenting a consistent interface.

### Usage

The relevant MPI functions are accessible through wrapper functions with
the underscore in the name removed. For example, the MPI routine
`MPI_AllReduce` is available as `MPIAllReduce`.

The arguments to the wrapper routines are a subset of those specified in
the MPI standard, and the explanations of those values may be found
there.

Many of the routines have versions which require specifying the array
dimensions manually (`*_len`) or automatically (`*_auto`). The latter
obtain the dimensions of the arrays to communicate by analysis of the
array bounds and are preferred, as they are less prone to user error.

If communication is desired between subsets of the available processors
(as is used in the CCMC code), the MPI communicator may be specified
through the `Node` parameter. No further detail is provided here about
how to construct these objects.

If the necessary MPI routine is not present in the `Parallel_neci`
module it will be necessary to add it. This can be awkward as discussed
below.

@warning
The `inplace` MPI functionality is present (although disabled) in NECI.
For the time being it should not be used due to serious bugs in the
interaction between ifort and OpenMPI.
@endwarning

@warning
Different implementations of MPI are not interchangeable at runtime,
even if they initially appear to be. For instance, code compiled using
OpenMPI on the ifort compiler will run on MPICH built with the gnu
toolset, but no communication will occur, and `nProcessors` independent
(identical) simulations will run.
@endwarning

@warning
In keeping with the notion of failing early, and failing loudly, outside
of the initialisation and cleanup code, most of the templated MPI
routines do not quietly report errors in status variables, but will kill
the entire simulation with a runtime error.
@endwarning

### Important variables

A number of important control variables are available for use within
NECI.

-   `nProcessors`<br>
    The total number of processors initialised in the MPI calculation.

-   `iProcIndex`<br>
    The (zero based) index of the current processor. This will be
    between 0 and `nProcessors-1`.

-   `root`<br>
    The (zero based) index of the root (head) processor. This should be
    tested in code using: `if (iProcIndex == root)`.

-   `nNodes`<br>
    The number of nodes initialised in the MPI calculation. In most
    cases this will equal `nProcessors`. These values will only differ
    if the MPI space is being subdivided into smaller nodes.

-   `iNodeIndex`<br>
    This (zero based) index gives the position of the processor in the
    current node. For most calculations this will be zero on all
    processors.

-   `bNodeRoot`<br>
    This specifies if the current processor is the root processor of a
    node. I.e. that this processor is responsible for macroscopic
    communication. For most simulations this is `.true.` on all
    processors.

### Implementation

The current implementation of the `Parallel_neci` module works extremely
well. It may, however, be necessary to modify it, either to deal with
changes in the available compilers, or (more likely) to add support for
additional MPI functions.

Fortunately, the latter is much more straightforward than the former! A
modified version of a currently implemented routine is likely to be the
best approach.

## Shared memory

In principle, when using MPI all of the processes are entirely
independent except for the explicit communication made through the MPI
library. Within NECI, as all of the particles are entirely independent
between annihilation steps this is a good thing.

However, we have a large amount of read-only data. The integrals which
are read in from the FCIDUMP files can consume a non-trivial proportion
of the available system memory, and are duplicated on each of the
processes. On modern multi-processor and multi-core systems this is an
outrageous waste of system memory — which is made worse by the fact that
if the same regions of this memory are requested on multiple different
processors the computer will not know these are the same, and will not
be able to make use of the L1, L2 and L3 cache to speed memory access.

The shared memory provisions in NECI substantially abuse the MPI
specification by mixing the available memory address space between the
different processes. This is extremely useful when done careful, but
there are some caveats to be aware of.

### How to use shared memory provisions

The shared memory provisions have been designed to be as interchangeable
with the normal Fortran memory allocation provisions. Essentially memory
is allocated largely as usual, but all of the processes on the same
physical node will end up with pointers to the same block of memory.

Memory is allocated with a call to

```Fortran
subroutine shared_allocate(name, ptr, dims)
```

This is a templated routine, and so will work for a wide variety of data
types and array dimensions. The `name` parameter specifies a unique name
for this array. This name is used to identify the particular block of
shared memory between the processes (there may be an arbitrary number).
The `ptr` parameter is a Fortran `pointer` with the relevant data type
and the correct array shape. `dims` is an array of integers indicating
the size of each dimension of the array to be allocated.

As an example,

```Fortran
real(dp), pointer :: real_arr(:,:)
...
call shared_allocate("example", real_arr, [6, 13])
```

will allocate a 2-dimensional array named “example” with the array
bounds `(1:6, 1:13)`, equivalent to `allocate(real_arr(6, 13))`.

The read only data should now be written. In principle, the data only
needs to be written from one of the processors per node — in practice it
is easiest to get all of the nodes to write the read only data. It is
important that data is only written at this stage — any processing that
needs to be done on this data must not happen in the shared memory, as
it could be disrupted by the other processes. This scheme only works
because all of the processes write exactly the same data.

After data initialisation, the processes must be synchronised with a
call to `MPIBarrier`. This will ensure that the read-only data is in the
same state in all of the threads.

@warning
The POSIX memory model makes few guarantees about *when* memory that is
shared between processes gets updated (as opposed to between threads).
This shared memory ***must not*** be used for communicating between the
MPI processes. That is what MPI is for.
@endwarning

@warning
There are no synchronisation primitives provided for use with these
shared memory blocks. No assumptions can be made about the order of
access between the processes. All actions that are carried out
***must*** be inherently thread-safe, or errors will eventually (and
unexpectedly) occur.
@endwarning

The data should be deallocated using the routine `shared_deallocate`.

### Quirks and limitations

-   **Customising data types**<br>
    The shared memory routines are templated so that they can be used
    with a wide variety of plain-old-data types and array sizes. For any
    data types outside of those already templated, new configuration
    options will be needed in `lib/allocate_shared.F90.template`. The
    information required is the number of bytes per element of the
    array.

    If custom `types` are desired, this may be quite tricky. It is
    difficult to guarantee the size, memory alignment and packing of the
    data in a custom type — the compiler is free to choose how it wishes
    to do this, and as such it varies from compiler to compiler. It is
    necessary to use the `sequence` keyword must be used to force the
    compiler to store data contiguously. If there are `pointer` or
    `allocatable` elements in the custom type, it cannot be used with
    shared memory, as the size of these elements is highly variable
    between compilers.

-   **Storing bit representations**<br>
    There is a special allocate routine `shared_allocate_iluts` which
    can be used for storing bit representations — these differ in that
    the lower bound of the first array index must be \(0\), rather than
    \(1\).

-   **Non-uniform memory architectures**<br>
    Modern chip architectures (in particular the intel i3, i5 and i7
    series of processors) move the memory controller onto the processor.
    This dramatically improves memory access speeds.

    The cost is that on multi-processor (as opposed to multi-core)
    architectures the memory is segmented into regions that are
    controlled by the different processors. Memory access within the
    region controlled by the current processor is faster than that
    between them.

    A performance increase could be obtained by only sharing memory
    between processes which are hosted on the same physical processor.
    This will require using the processor affinity options of MPI, which
    is not routinely done currently with NECI. Further, it will require
    subdividing the processes on each physical node into sub-nodes —
    support for this exists in the code in NECI, but the processor
    specific information required to correctly subdivide the processes
    is not included, and this facility is currently switched off.

    For processor topology, see the Intel documentation at
    [https://software.intel.com/en-us/articles/
                intel-64-architecture-processor- topology-enumeration](https://software.intel.com/en-us/articles/
                intel-64-architecture-processor- topology-enumeration).

-   **Disabling shared memory**<br>
    Shared memory is enabled by the pre-processor define `SHARED_MEM_`.
    This can be found in the config files for platforms that support it.
    If necessary, this can be removed from relevant config files and
    Makefiles, which will disable inter process memory sharing, and fall
    back on the normal Fortran `allocate` mechanism.

### How it works

The way that inter-process memory sharing works is system specific, and
depends on operating system primitives being used directly.

The code templating functionality in NECI is used to make the shared
memory wrapper work for a wide range of data types, with differing
numbers of dimensions in the arguments. From this the size of memory
required in bytes is calculated, and this is passed to a C helper
function that allocates the memory. The returned raw pointer is
converted to a Fortran one through the Fortran 2003 intrinsic
`c_f_pointer`.

The C helper library contains a number of utility functions to help it
interface with Fortran. In particular wrapper so it can print to the
same output, and a mapping to store additional data about the allocated
memory to assist in deallocating without the Fortan code requiring
knowledge of how the operating system intrinsics work.

Beyond that, there are three main code paths:

-   **POSIX**<br>
    A unique file name is created, based on the current working
    directory. The function call `shm_open` with the control parameter
    `O_CREAT` will open, or create, a POSIX shared memory object. As a
    result, one of the processes on a node will create the mapping, and
    the others will join it — it is unknown in advance which will do the
    creating.

    The returned descriptor acts as a file, and so is set to the desired
    length with `fdtruncate`, and then mapped into memory as if it were
    a memory-mapped file using `mmap`. The file descriptor is then
    closed as it is no longer needed.

    Once all of the processes have mapped the memory, the shared memory
    object is unlinked, ensuring that the memory will be deallocated by
    the operating system when all of the processes stop (preventing a
    memory leak if the processes crash).

    The memory can then be manually deallocated using `munmap`.

-   **System V**<br>
    A unique file name is generated, based on the current working
    directory, and the file is created. A System V Inter Process
    Communication Key is created and obtained from this file using the
    routine `ftok`.

    Using this key, a shared memory object of the correct size is
    created using `shmget`, and this region is mapped into memory using
    `shmat`.

    Once all of the processes have reached this point, the shared memory
    control object is destroyed using `shmctl` to ensure the operating
    system will deallocate the memory in case of a crash.

    The memory can be manually deallocated using `shmdt`.

-   **Windows**<br>
    The memory mapped file provisions in Windows are used to generate a
    shared memory region. A unique file name is created, based on the
    current working directory, and then a non-filesystem backed file is
    created with `CreateFileMappingW`. This can be opened by all of the
    other processors on the same node using `OpenFileMappingW`, and this
    "file" is mapped into memory on all of the processes using
    `MapViewOfFile`.

    The memory can be manually deallocated using `UnmapViewOfFile`
    followed by `CloseHandle`.

The subroutine `iluts_pointer_jig`, which is used when allocating bit
representations which start from an index of \(0\) rather than \(1\), is
an interesting demonstration of how to manipulate the bounds of arrays
declared in Fortran in compilers that do not support the array reshaping
assignments described in Fortran 2003 (most compilers).

## Integral retrieval

Integral retrieval is found in the tightest of the tightest loops within
NECI. Calculating each Hamiltonian matrix element may require a number
of different two- and four-index integrals, and at least one Hamiltonian
matrix element is required for each generated excitation.

As a result, the normal approach taken to prevent code repetition
introduces a bottleneck. If there is a `get_umat_el` function, this will
have to contain the logic as to which type of integral is being
obtained, and where this should be located. This conditional logic will
be executed for every required integral.

As such, access to the integrals is through a procedure pointer, with
the interface

```Fortran
abstract interface
    function get_umat_el(i, j, k, l) result(hel)
        use constants
        implicit none
        integer, intent(in) :: i, j, k, l
        HElement_t :: hel
    end function
end interface
```

This function pointer can then be called from anywhere in the code
(after it is initialised) and this will call the correct routine.

### FCIDUMP files

The most commonly usedintegrals routines store their integrals in memory
after reading them in from a `FCIDUMP` file.

These routines support FCIDUMP files produced by various codes,
including MOLPRO, PSI3 and VASP (hacked versions of Dalton and QChem
also support this - I believe MOLCAS as well now). Certain (generally
legacy) FCIDUMP files have a strictly fixed format, which can result in
adjacent columns of indices merging. As such, they are read via an
explicit format. To use MOLPRO or other FCIDUMP files the `MOLPROMIMIC`
option should be supplied in the `SYSTEM` block of the input file. Other
sources of FCIDUMP file should use the `FREEFORMAT` option. In general,
this should *always* be used.

The FCIDUMP files can be briefly summarised by

-   **Header**<br>
    The header is specified as a Fortran namelist, with the following
    possible elements:

    -   `NORB`<br>
        Specifies the number of spatial orbitals in the system if RHF,
        or the number of spin-orbitals if UHF.

    -   `NELEC`<br>
        Specifies the number of electrons the Hartree–Fock determinant
        should have.

    -   `MS2`<br>
        Specifies twice the total projected spin of the system (such
        that it is always an integer)

    -   `ORBSYM`<br>
        A list of spatial symmetries of the orbitals specified in (???
        Check name with Giovanni) format in orbital order. Note that the
        structure of the reference determinant will be visible here if
        the `MOLPROMIMIC` option is used, and the reference is to be
        determined by the order of the orbitals.

    -   `ISYM`<br>
        The symmetry of the reference determinant specified in the same
        format.

    -   `UHF`<br>
        T if the file contains a UHF basis, otherwise F.

    -   `SYML, SYMLZ`<br>
        The total and projected orbital angular momentum for systems
        with an axis of rotation. Currently FCIDUMP files which use this
        option can only be generated using QChem, and the resultant file
        must be pre-processed using the TransLz utility.

    -   `PROPBITLEN, NPROP`<br>
        These parameters are used to describe the behaviour of k-point
        symmetries. For more details contact George Booth.

-   **4-index integrals**<br>
    Following the header, all of the integral and energy lines may
    follow in arbitrary order.

    The 4-index integrals are specified with the format `Z i j k l`,
    where `Z` is the `real` or `complex` element, and the indices are
    integers. All indices are one-based. Indices are in **chemical
    notation**!

-   **2-index integrals**<br>
    The 2-index integrals are specified in the same way with the final
    two indices equal to zero.

-   **Fock energies**<br>
    The Fock energies are specified in the same way, with the final
    three indices equal to zero.

    These values are used for determining the initial (Hartree–Fock)
    determinant, which is constructed from the lowest energy
    spin-orbitals. They are strictly optional. If the `MOLPROMIMIC`
    option is supplied then the order of the orbitals determines the
    reference determinant.

-   **Core energy**<br>
    The core energy is specified in the same way, but with all four
    indices set to zero.

### Generation on the fly

A number of schemes exist for generating integrals on the fly. In
particular the Uniform Electron Gas and Hubbard model have integrals
that can be trivially calculated, and so FCIDUMP files are not used.

### Nested schemes

There exist two function pointers with the same abstract interface,
`get_umat_el` and `get_umat_el_secondary`. Once the primary method of
determining integrals has been selected, then this be moved to the
secondary pointer, and another wrapper routine used instead. This allows
additional filtering logic to be applied.

The wrapper routine should then call `get_umat_el_secondary` internally.

### Fixed Lz and complex orbitals

If either projected angular momentum (Lz) symmetry, or complex orbitals
are being used, the pattern of restricted zeros in the integrals
changes. Wrappers around the normal `get_umat_el` routine are used to
supply these zeros.

These are examples of the nested schemes above.

### Further schemes

It is likely that further schemes will be required in the future. In
particular, approximate, and interpolating schemes are likely to be
required to reduce the \(\mathcal{O}(M^4)\) memory dependence on the
number of orbitals. These should be implemented as new routines to be
pointed at using the function pointers.

## Hamiltonian matrix element evaluation

The two- and four-index integrals are the primary input to a FCIQMC
calculation. However, they are used via Hamiltonian matrix elements.
There are a number of considerations for the way that Hamiltonian matrix
elements are used in NECI.

### Different ways to obtain Matrix elements

Certain pieces of information are required in order to calculate
Hamiltonian matrix elements. The Excitation level (the number of
differing orbitals between the determinants), the corresponding
excitation matrix and the parity of the excitation are necessary.
Depending on the excitation level, the decoded list of occupied orbitals
is required.

Depending on where in the execution path the Hamiltonian matrix elements
are required, different information is readily availabile. As some of
this is relatively expensive to calculate (in particular the parity) it
is important that as much information as is available is used.

As part of the spawning and particle generation process, there is a
function pointer called `get_spawn_helement`. This is passed the natural
integer and bit representations of the determinant, the excitation
level, excitation matrix, parity and (potentially) pre-calculated matrix
element. Depending on the initialisation options, the routine which is
selected will make use of the correct subset of these. It is important
that the excitation generator and the choice of Hamiltonian matrix
element generation here are tightly coupled.

Elsewhere in the code, the routine `get_helement` is used to obtain
matrix elements.[^gethelement] This routine comes in a number of flavours. The
available options are

-   `nI, nJ`<br>
    This routine will return the matrix element between any two,
    arbitrary, decoded determinants.

-   `nI, nJ, ic`<br>
    This version is provided the excitation level of the determinant in
    addition.

-   `nI, nJ, ic, ilutI, ilutJ`<br>
    The bit representations of the two determinants are provided, which
    greatly enhances calculating the parity if needed.

-   `nI, nJ, ilutI, ilutJ`<br>
    If both the bit representations and the decoded versions are present
    but no further information is known then this form should be used.

-   `nI, nJ, ilutI, ilutJ, ic_ret`<br>
    This version is the same as the above, but returns the excitation
    level of the pair of determinants in addition to the matrix element.

-   `nI, nJ, ic, excitMat, tParity`<br>
    When everything about the relationship between the two determinants
    is fully known, this form should be used. It is used implicitly
    after excitation generation when the excitation level, matrix and
    parity are known. The second decoded determinant is not used in
    determinental calculations, but is provided to support systems such
    as CSFs.

For calculating diagonal Hamiltonian matrix elements, the routine
`get_diag_helement` should be used[^diagelem]

#### The cost of the parity

The overall sign of the returned Hamiltonian matrix element is modulated
by the parity of the excitation — that is, the whether the number of
pairwise swaps of orbitals required to maximimally align the two
determinants according to the standard order of the first is odd or
even.[^double_excit]

The parity may be obtained by actively aligning the decoded
representations, or by examination of the bit representations. The
latter is much faster than the former, but is still the rate limiting
factor for calculating Hamiltonian matrix elements.

Where this factor is known (i.e. after excitation generation) it should
be used. It is worth noting that the weighted excitation generation
scheme is only viable because only the absolute value of the matrix
elements is modelled, so the parity generation step can be entirely
ignored.

#### Slater–Condon rules implementations

The Slator–Condon rules are implemented in the file `sltcnd.F90`, in the
routines `sltcnd` and then more specifically `sltcnd_0`, `sltcnd_1` and
`sltcnd_2`. There is a certain amount of duplication both of code paths
and of conditional testing in these routines — which is a clear
violation of the DRY principle.

Access to matrix elements is inside the smallest of the tight loops in
NECI. The provision of duplicated logical pathways, rather than
conditional switching, has a measurable performance benefit in this
case.

#### HPHF matrix elements

In most cases the HPHF matrix elements are trivially obtained from the
normal determinental matrix elements through multiplication by a factor
of \(\sqrt{2}\) in all cases where the determinant is not closed shell.

If the excitation level is less than or equal to 2, then the elements
are a little more complicated. It is possible that both the target
determinant and its spin pair are connected to the source determinant,
and so the resultant matrix elements will require two calls to the
Slater–Condon routines.

#### Spin eigenfunctions

Spin eigenfunctions may be expressed as linear combinations of all
determinants with a given spatial structure. The Hamiltonian matrix
elements may be expressed as a two index sum over all the determinants
with each of the involved spatial configurations.

Given that the number of configurations increases permutationally with
the number of unpaired electrons, this scheme scales extremely badly.
NECI has some aggressive optimisation to slightly improve this scaling
(see Simon’s thesis), but ultimately it is a dead end for big
calculations.

See the SPINS branch for Hamiltonian matrix element calculation between
other types of spin eigenfunctions. These can scale extremely well, but
introduce other problems.

## Excitation generation

Excitation generation is the most complicated and intricate part of
NECI. Not only must random moves be made, but they must be made in a
manner which is efficient (taking account of symmetry, and in terms of
implementation). They must also correctly compute the generation
probabilities, taking into account multiple possible selection routes to
the same determinant.

A failure to generate all of the connected determinants, or mistakes
made in calculating the generation probabilites will lead to silent
errors that manifest themselves only through (potentially subtly)
incorrect energies being produced by the simulation.

There are a large number of factors that must be considered when writing
an excitation generator. This section aims to give an overview of some
of them, and a brief outline of the two most commonly used excitation
generators.

The Alavi group collectively has a large amount of experience writing
excitation generators, and anybody intending to write or modify them
would be advised to speak to Simon Smart or George booth first.

### Interface

Excitation generation is accessed through a procedure pointer. As such,
all excitation generators *must* have the same signature. It is
described by

```Fortran
subroutine generate_excitation_t (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                  ex, tParity, pGen, hel, store)

    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use FciMCData, only: excit_gen_store_type
    use constants
    implicit none

    integer, intent(in) :: nI(nel), exFlag
    integer(n_int), intent(in) :: ilutI(0:NIfTot)
    integer, intent(out) :: nJ(nel), ic, ex(2,2)
    integer(n_int), intent(out) :: ilutJ(0:NifTot)
    real(dp), intent(out) :: pGen
    logical, intent(out) :: tParity
    HElement_t, intent(out) :: hel
    type(excit_gen_store_type), intent(inout), target :: store

end subroutine
```

`nI` and `ilutI` describe the natural integer and bit representations of
the source determinant for the excitation. `nJ` and `ilutJ` will store
the generated excitation. `exFlag` specifies the type of excitation —
this is generally unused, but in some excitiaton generators permits
selecting between single and double excitations for the purposes of
testing.

`store` provides a location for the excitation to store information that
is specific to the source determinant. This information is available for
further excitations from the same site.

`ex` will contain the excitation matrix; `ex(1,:)` contains the
*spin-orbitals* of the source electrons in `nI`, and `ex(2,:)` contains
the *spin-orbitals* that have replaced these values.

`ic` is used to return the excitation level of the resultant site
relative to the source site, `pGen` returns the generation probability
for the generated determinant, `tParity` returns the parity of the
excitation (is the number of pairwise swaps required to align `nI` and
`nJ` even or odd), and `hel` can return the Hamiltonian matrix element
for this excitation so that it need not be calculated later (this is
only really of interest for HPHF calculations, where it can be easier to
calculate this matrix element in the excitation generator for technical
reasons).

@warning
Note that the determinant returned in `nJ` *must* be sorted in the same
fashion as the source determinant `nI`. This normally means in order of
ascending spin-orbital number. This can result in substantial
reshuffling of the internal order of the detrminant.
@endwarning

@warning
It is possible that to implement a sensibly efficient excitation
generator, the simulation may make choices which leave no possible
excitations remaining. An excitation may be aborted by setting the first
element of `nJ` to zero, so long as the calculated generation
probabilites for successful excitations are correct.
@endwarning

### Symmetry handling

The excitation generator must always preserve the symmetry of the
determinant. The way in which this is done depends on the excitation
generator, but will involve measuring the symmetry of the source
orbitals and selecting the target orbitals appropriately.

The symmetry of a specific orbital, `orb`, may be found in the array
`G1`. The spatial symmetry is found at `G1(Orb)%Sym%S`, the k-vector (if
appropriate) at `G1(orb)%k`, and the projected spin and orbital angular
momentum values at `G1(orb)%Ms` and `G1(orb)%Ml`.

There are also a number of useful macros for measuring and manipulating
spin values. `is_beta` and `is_alpha` test if orbitals have beta and
alpha spin respectively, whilst `is_one_alpha_beta` test that two
orbitals have differing spins. The `get_spin` macro gets the spin in the
unusual format for \(1\) for alpha and \(2\) for beta as used in some
other NECI routines, and `get_spin_pn` obtains the more standard values
of \(\pm 1\). The `is_in_pair` macro tests if two orbitals belong to the
same spatial orbital (including the case where they are the same
orbital).

The macros `get_alpha` and `get_beta` obtain the alpha or beta orbital
corresponding to the same spatial orbital as the current spin orbital
(which may or may not be the same orbital). `ab_pair` obtains the spin
paired orbital to the current one.

It is normally necessary to combine and manipulate symmetries, as for a
double excitation \(\Gamma_i\otimes\Gamma_j = \Gamma_a\otimes\Gamma_b\),
but there is no reason for any of \(Gamma_i,Gamma_j,Gamma_a,Gamma_b\) to
be the same. The products of symmetry labels should be taken using
`RandExcitSymLabelProd`, rather than directly using `ieor` commands
(which will normally obtain the same result), to protect against the
consequences of using different types of symmetry. When combining Ml
values it is important to bear in mind the maximum permitted value of
Ml, and abort the excitation if necessary.

The class count arrays, described above, are indexed by a a combined
symmetry index, that combines spin, k-point, momentum and spatial
symmetries. The index to this array is provided by the function
`ClassCountInd`. Given a total spin, momentum and spatial symmetry, the
paired class count index may be obtained using `get_paired_cc_ind`. The
number of occupied, and unoccupied spin orbitals are stored in the data
store as decribed above. The total number of orbitals corresponding with
a given class count is found at the corresponding location in the array
`OrbClassCount`.

The array `SymLabelList2` contains all orbitals sorted by symmetry
class, with offsets given in `SymLabelCounts2`, such that all of the
orbitals with the same symmetry as a given orbital, `orb`, may be found
as follows

```Fortran
integer :: sym, spn, ml, norb, cc_ind, offset

! Obtain the symmetry labels associated with this orbital
sym = G1(orb)%Sym%s
spn = get_spin(orb)
Ml = G1(orb)%Ml

! Obtain the combined symmetry index
! cc_ind = ClassCountInd(orb) ! ... alternatively
cc_ind = ClassCountInd(spn, sym, ml)

! Get position and count of orbitals
norb = OrbClassCount(cc_ind)
offset = SymLabelCounts2(1, cc_ind)

! And this slice contains the appropriate orbitals (including orb)
SymLabelList2(offset : offset + norb - 1)
```

### Manipulating determinants, and bit representations

The primary output of the excitation generator is the newly generated
determinant. It is much more computationally efficient to copy and
modify the source determinant and bit representation than to start from
scratch. As such there are two processes to consider.

-   **Manipulating bit representations**<br>
    Orbitals can be cleared from the bit representation using the macro
    `clr_orb`, and set using the macro `set_orb`. For an excitation from
    orbital `a` to orbital `i` this would look like
    ```Fortran
    clr_orb(ilut, a)
    set_orb(ilut, i)
    ```

-   `make_single and make_double`<br>
    Manipulating the natural integer representation of determinants is
    substantially more complicated, as there is a *requirement* that
    they remain sorted by spin-orbital number. Although it would be
    straightforward to, essentially, replace the excited orbital with a
    new one, and sort it, this is inefficient as the parity of the
    excitation is going to be required. If the manipulation can be
    carried out in such a way that the parity is naturally returned this
    will substantially improve the efficiency.

    The routines `make_single` are passed the source determinant and
    respectively one or two source electron positions and target
    orbitals. They perform the substitution and measure how much the
    newly placed orbitals must be moved by to restore sorting
    (accounting for edge cases).

    These functions return the sorted determinant, the excitation matrix
    and the parity of the excitation. They should be used by all
    determinental excitation generators.

### Data storage

Each site in the main particle list may be occupied by multiple
particles. There is a reasonable amount of information which the
excitation generator must generate for each site it excites from — it
makes sense to store this information and reuse it for each of the
particles on the site. A data structure of type `excit_gen_store_type`
is passed to the excitation generator that it can use for this purpose:

```Fortran
type excit_gen_store_type
    integer, pointer :: ClassCountOcc(:) => null()
    integer, pointer :: ClassCountUnocc(:) => null()
    integer, pointer :: scratch3(:) => null()
    integer, pointer :: occ_list(:,:) => null()
    integer, pointer :: virt_list(:,:) => null()
    logical :: tFilled
    integer, pointer :: dorder_i (:) => null()
    integer, pointer :: dorder_j (:) => null()
    integer :: nopen
end type
```

The arrays that are required for the current excitation generator are
allocated during calculation initialisation by the routine
`init_excit_gen_store`.

When the excitation generator is being called on a new determinant, the
`tFilled` variable will be set to `.false.` by the main loop code. Once
the generator has filled in the required information this should be set
to `.true.` to indicate that the data may be reused.

The arrays `ClassCountOcc` and `ClassCountUnocc` are initialised by
calling `construct_class_counts` at the start of the excitation
generator. They contain a count of the number of occupied, and the
number of unoccupied, orbitals associated with each spin-symmetry
combination. These are used for calculating the probabilities in
essentially all excitation generators.

The `occ_list` and `virt_list` variables store only information relevant
to the excitation generator in `symrandexcit3.F90`.

The `dorder_*`, `nopen` and `scratch3` variables are used in the CSF
excitation generation code.

As there is only one instantiation of this data structure, in the main
loop, adding additional elements adds only negligible overhead. If
future excitation generators need a place to store information between
calls, it should be added here.

The excitation generation routines should be as tightly coupled to the
Hamiltonian matrix element generation routines as possible. In
particular, information such as the excitation level and the excitation
matrix are trivially available in the excitation generator and expensive
to calculate otherwise. As it stands, the interface for the excitation
generator permits passing the excitation level, parity and excitation
matrix to the matrix element generation routines — this is sufficient
for determinental systems.

If NECI is to be extended to efficiently consider other basis functions,
then other data will be required. A structure should be created to pass
this information around, so that each element which is added does not
require adjusting the interface of each and every excitation generator.
This has been done in the SPINS branch, with a type named
`spawned_info_t`, which is unlikely to be merged into master for
independent reasons, but may provide a template for doing this if
required in the future.

### Re-use and rescaling of random numbers

In determinental calculations, even with large basis sets, each site is
connected to at most a few thousand others. A single 64-bit random
number contains vastly more random information than is required to make
a good choice between all of these options.

However, most of the potential algorithms for excitation generation
involve a hierarchy of choices. The natural implementation of these
choices is to generate a new random number for each new choice that is
required. This is highly wasteful, especially as random number
generation is relatively expensive.

More bang-for-your-buck can be extracted from random numbers by
partitioning them, and rescaling them as you move through the hierarchy
of choices. A simple example demonstrates this, considering the
macroscopic choice between single excitations and double excitations:

```Fortran
real(dp) :: r
r = genrand_real2_dSFMT()
if (r < pDoubles) then
    ! Double excitation selected. Now rescale r
    r = r / pDoubles
    ...
else
    ! Single excitation selected. Now rescale r
    r = (r - pDoubles) / (1.0_dp - pDoubles)
    ...
end if
```


This is not always done in the current implementations of excitation
generators, but should always be considered for efficiency.

@warning
This technique should not be used where selections may need to be
redrawn. Repeated repartitioning and scaling of the random number will
rapidly deplete the available randomness if it is repeated a number of
times due to redrawing.
@endwarning

### Timestep selection and other control parameters

The spawning process is normally the limiting factor for the selection
of the imaginary timestep. The magnitude of the spawn, \(n_s\), is given
by
\[n_s = \delta\tau \ensuremath{\left\lvert\frac{H_{ij}}{p_\mathrm{gen}(j|i)}\right\lvert}.\]
As such, a limit on the value of \(\delta\tau\) may be obtained from the
maximum value of that ratio, combined with a (chosen) maximum size of
spawn;
\[\delta\tau_{max} = n_{s,max} \times \left(\max\ensuremath{\left\lvert\frac{H_{ij}}{p_\mathrm{gen}(j|i)}\right\lvert} \right)^{-1}.\]
The ratio should be accumulated through a calculation, and can be used
to determine the optimum time step on the fly.

This approach may be generalised to a larger number of parameters.
Additional control parameters that influence decisions in the excitation
generator may be introduces (such as the choice between single and
double excitations, or the extent to which a bias is made towards or
against double excitations which are spin aligned).

Considering the different choices that are made in the excitation
generator, this splits the excitations into a number of different
categories. For optimal timestep values the worst case for each of these
categories should be optimised to give the maximum spawn size.

The impact of each of the parameters to be optimised must be removed
from the generation probability values, so that an unaffected value can
be maximised. For example considering a split between single and double
excitations; \[n_{s,max}
        = \delta\tau \frac{\delta\tau}{p_{single}}
            \ensuremath{\left\lvert\frac{H_{ij}p_{single,iter}}{p_\mathrm{gen}(j|i)}\right\lvert}_{single}
        = \delta\tau \frac{\delta\tau}{1.0 - p_{single}}
            \ensuremath{\left\lvert\frac{H_{ij}(1.0 - p_{single,iter})}{p_\mathrm{gen}(j|i)}\right\lvert}
            _{double}.\] In this case the current value of
\(p_{single,iter}\) (i.e. the value on the particular iteration the
generation occurred) is removed from the ration which can then be
maximised. Closed form expressions for the optimal values of
\(\delta\tau\) and \(p_{single}\) can then be found straightforwardly.

This approach can be extended, although the means to isolate the impact
of parameters on the generation probability, and the structure of the
final closed form expressions will depend on the form of the excitation
generator.

### Testing excitation generators

It is critical that excitation generators are *correct*. That is that
all of the connected determinants are generated, and that the generated
probability is correct.

There are a number of metrics that can be used to test this. A testing
function should be written that takes a given source determinant as a
parameter, and runs the excitation generator a very large number of
times on this determinant (for example, 10 million times). This test
function should contain a number of tests:

-   **Are the correct determinants generated**<br>
    A list of all single and double excitations of the correct symmetry
    should be enumerated in a brute force manner (see
    `GenExcitations3`). It should then be checked that *all* of the
    determinants in this list are generated by the excitation generator,
    and no determinants outside this list are generated.

-   **Normalisation of generation probabilities**<br>
    An accumulator value should be kept for each possible determinant
    that can be generated. Each time the determinant is generated, the
    value \(p_\mathrm{gen}^{-1}\) should be added to that determinants
    accumulator.

    On average this will add a value of \(1.0\) to the accumulator for
    every excitation generation attempt. Thus, when all of the
    accumulated values are divided by the number of generation attempts
    made by the test function, all of the values should give roughly
    \(1.0\).

    There may be a reasonable amount of variation, but repeated across a
    few different runs with different random number seeds it should be
    clear that all of these values give roughly \(1.0\).

    Note that if some connections are *extremely* strongly weighted
    against (as can happen with the weighted excitation generator) they
    will sum in a very large term very rarely, and so their averaged
    value can jump around quite dramatically.

-   **Overall probability normalisation**<br>
    As an extension to the above, the sum of all of the accumulators
    should stocastically tend towards the number of connections
    multiplied by the number of spawning attempts.

    If this total accumulated value is divided by the number of
    connections multiplied by the number of spawning attempts it should
    average very strongly to \(1.0\), normally to four or five decimal
    places.

For an example test function, see `test_excit_gen_4ind` in
`symrandexcit4.F90`. An excitation generator that passes these tests is
not guaranteed to be correct, but it is quite likely.

### How this interacts with HPHF functions

HPHF functions present a complexity for excitation generation.
Excitation generation occurs on determinants, but the HPHF function
includes the spin flipped pair. When there are only a small number of
unpaired electrons, the possibility of the excitation generator having
generated the spin flipped version also needs to be available, although
this is not readily accessible in the excitation generator.

As a result, a routine must be provided to calculate the probability of
generating \(\ket{D_j}\) from \(\ket{D_i}\) without actually generating
an excitation.

@warning
Once the excitation generator has been planned out, it makes sense to
write this calculation function first (it tends to be simpler than the
excitation generator). A call to it can then be included at the end of
the excitation generator inside `#ifdef DEBUG_` flags, which provides a
powerful correctness check, and catches any cases where these functions
diverge.
@endwarning

### The uniform selection excitation generator

A full discussion of this excitation generator is found in the paper
(...). The operation of this excitation generator is also the basis for
a number of other more specialised excitation generators (such as the
CSF excitation generator). This excitation generator is implemented in
`symrandexcit2.F90`.

Although the generation probabilities associated with this excitation
generator are highly non-uniform, it fundamentally makes choices at each
*stage* of the excitation generator uniformly between the options
presented at this point.

For single excitations, an electron is chosen uniformly at random. This
fully specifies a symmetry index, and from the appropriate list a vacant
orbital is chosen at random. This selection is generally made by picking
an orbital with the appropriate symmetry at random, and re-drawing if it
is already occupied (unless there are very few available choices, in
which case the \(n\)-th available orbital is chosen from an enumeration
of available orbitals).

For double excitations, a pair of electrons are chosen at uniform (using
a triangular mapping). This specifies the total symmetry of the
excitation. A vacant orbital, \(a\), is then chosen at random from the
entire array of orbitals (in the same way as above, redrawing if an
occupied one is chosen, or if an orbital is selected that has no
available orbital \(b\) that could produce the correct total symmetry).
This in turn specifies the symmetry of orbital \(b\), which is chosen in
the same way as the orbital for a single excitation.

The possibility of having picked the orbital \(b\) first, and then \(a\)
must be accounted for in calculating the generation probabilities.

### The weighted excitation generator

This excitation generation scheme is discussed in a great deal more
detail in a paper which is currently pre-publication. This excitation
generator is implemented in `symrandexcit4.F90` as
`gen_excit_4ind_weighted`.

Fundamentally, a weighting value is used to bias the choice of orbitals
towards those which are likely to correspond to large Hamiltonian matrix
elements. This is done using a Cauchy–Schwarz decomposition of the
Hamiltonian matrix elements. The specific elements used depend on the
spin choices made, and are discussed in the referred paper.

For speed, the entire choice of orbitals of a given symmetry are
considered, and the weightings corresponding to orbitals which are
occupied in the source determinant are set to have a zero weighting.

For single excitations, an electron is chosen uniformly at random. This
determines the target symmetry. A cumulative list of the weighted terms
for each of the available orbitals is generated, and a random number on
the range of the largest element in this list is generated. A specific
index into this list is selected by binary searching for the first
element in the cumulative list which is greater than or equal to this
random number, and this index specifies which of the orbitals of the
given symmetry are selected as the target orbital.

For double excitations, a choice is made of whether the two orbitals
should have the same spin or different spin (the probability of either
is optimised on the fly at run time). If the two electrons are to have
the same spin, they are chosen uniformly using a triangular mapping, and
otherwise uniformly using a rectangular mapping. This entirely
determines the target combined symmetry index.

A list is constructed, with a length equal to the total number of
available target symmetry indices. This is populated with a cumulative
count of the product of the number of available orbitals of a given
symmetry with the number of available orbitals of the paired combined
symmetry required to result in the determined total symmetry. To prevent
double counting, there are considered to be no combinations where the
first index is larger than the second index. A pair of symmetries are
chosen using binary searching as described before, weighted towards the
symmetries with the most available choices.

Two orbitals, \(a,b\), with the corresponding symmetries, are then
picked in the same way as the orbital was selected for single
excitations (although the weighting terms must now consider the effects
of both source orbitals). When choosing orbital \(b\), if it comes from
the same symmetry category as orbital \(a\) the relevant elements must
be adjusted in the cumulative list.

For calculating the probabilities, the only duplication of generating
the same determinant is for the case when both orbitals have the same
spin and symmetry. This must be accounted for in the generation
probability.

## Determinant data storage

There are two main locations where blocks of information concerning
particles associated with determinants are stored.

### Transmitted data

This concerns data that will be (or could be) transmitted between
different processors. Essentially this includes persistent information
that is generated by the excitation generator. Combined, this
information forms the *bit representation* of a determinant.

All data access should be through the getter and setter functions in the
`bit_reps` module (`BitReps.F90`), with the exception of the orbital
representation which may be manipulated directly as it will always be
first. There are general `encode_bit_rep` and `extract_bit_rep` routines
which package and extract all of the relevant information, as well as
specific accessors.

In most of the code the encoded values for specific determinants are
referred to as `ilut`, standing for “Integer Look-Up Table” which
relates to the orbital representation.

The packaging of this data varies substantially between different
compiler and runtime configurations. For example, integer runs on 64-bit
machines co-opt some of the additional bits in the storage of the signed
particle counts to store the flags. It is essential that no attempt is
made to access the data in this array direcly.

The representation of each determinant is of dimensions `0:nIfTot`, that
is of length `nIfTot+1`. Unusually for an array in Fortran it is a
\(0\)-based index. The component of the overall representation required
to uniquely determine one site is `0:nIfDBO` (DBO stands for DetBitOps),
which includes the orbital description and any CSF descriptiors.

Each of the components of the representation has a `nIf*` (Number of
Integers For) and `nOff*` (offset of) value. These should not be used
directly unless adding data to this representation.

-   **Orbital description**<br>
    Determinants are stored by a bit-representation of the choice of
    orbitals to construct their Slater Determinants from. This
    representation is of length \(2M\), i.e. the same as the number of
    available spin-orbitals. An occupied orbital is represented with a
    \(1\), and a vacant orbital by a \(0\).

    This representation can be accessed directly. To simplify things,
    various macros are available. The macros `IsOcc` and `IsNotOcc` test
    the occupancy of particular orbitals, and the macros `set_orb` and
    `clr_orb` modify it without requiring the developer to worry about
    the data representation.

    To decode the bit representation of the orbitals into the natural
    orbital representation use the `decode_bit_det` routine.
    Equivalently, the `EncodeBitDet` generates the orbital bit
    representation from the natural integer version.

-   **CSF descriptors**<br>
    If CSF descriptor labels are required, they are stored here.

-   **Signed particle count**<br>
    The representation of the coefficient, or “sign” of the determinant
    is the primary value which is evolved during an FCIQMC simulation.

    These values are represented by an array of `real(dp)` coefficients
    of length `lenof_sign`. In a standard simulation this length is
    unity, but for complex walkers two values are used and the
    double-run code uses extras.

    These values can be used via the routines `extract_sign` and
    `encode_sign`. The specific elements of the sign array can be
    accessed by `extract_part_sign` and `encode_part_sign`. If all
    particles are to be removed the routines `nullify_ilut` and
    `nullify_ilut_part` can be used.

-   **Flags**<br>
    Any number of flags may be associated with a particular site. A
    current full list may be found in the file `bit_rep_data.F90`. The
    most commonly used flags are the `flag_is_initiator` and
    `flag_parent_initiator` flags — be aware that these are different
    names for the same thing (depends if considering particles in the
    main or spawned lists).

    There is a pseudo-flag, `flag_negative_sign`, that only exists to
    help combining the representation of the coefficients and flags. It
    should never be used directly outside of the representation
    manipulation routines.0

    These flags may be tested using the `test_flag` and `set_flag`
    routines referencing the specific flag desired from the above list.

    Further, then entire set of flags may be extracted or stored using
    the `extract_flags` and `encode_flags` routines. Once they have been
    extracted, they may be examined and manipulated with the `btest`,
    `ibset` and `ibclr` in the same way as the set and test functions
    above. A special routine `clear_all_flags` exists for resetting the
    flag status.

    It is extremely important that the flags element of the bit
    representation is not accessed directly.

### Locally stored data

This concerns information that will *never* be transmitted. All of this
information is either used in tracking the status of an occupied
determinant, or for values that could in principle be regenerated when
required but are stored for optimisation.

All data access should be through the `get_*` and `set_*` routines found
in the module `global_det_data`.

This data is stored in the global array `global_determinant_data`, but
this array should not be accessed directly! The location of data within
this array is determined by the initialisation routine in the module
`global_det_data`, and depends on the calculation being performed.

The data currently stored in this array are

The diagonal Hamiltonian matrix element for the determinant (`diagH`),

The average signed occupancy of the determinant (`av_sgn`, only for
RDMs).

The iteration on which determinants became occupied, (`iter_occ`, only
for RDMs) and

The moment in imaginary time on which the (contiguous) occupation of
this determinant began (`tm_occ`, only with experimental initiators).

## Transcorrelated integrals

The usage of the transcorrelated approach for ab-initio systems requires
us to include 3-body interactions, that come with 6-index integrals to
be handled by NECI, as well as triple excitations to be generated. The
triples excitation generator is uniform and is located in
`tc_three_body_excitgen.F90`, while the handling of the 6-index
integrals is done by the `LMat_class.F90` (for the storage and reading
of integrals) and the `LMat_mod.F90` (for getting matrix elements).

When reading 6-index integrals, they can either be read from an ASCII
formatted file, or an HDF5 file. When reading them from an ASCII file,
the file has to be formatted in a general FCIDUMP fashion, containing
all non-zero integrals \(L_{ijk}^{abc}\) in the format

```Fortran
<integral> i j k a b c
```

where `i`, `j`, `k` are the indices of the orbitals to excite from and
`a`, `b`, `c` those of the orbitals to excite to. No further information
is required in the file. The default filename is `TCDUMP`.

When reading the 6-index integrals from an HDF5 file, all data is
required to be contained in a group called `tcdump`, containing an
attribute named `nInts` containing the number of non-zero integrals, and
two datasets called `values` and `indices`, containing the values of the
non-zero 6-index integrals and their indices. The `indices` dataset has
to be of 6 times the size than the `values` dataset, each group of 6
indices is attributed to one value (in storage order).
