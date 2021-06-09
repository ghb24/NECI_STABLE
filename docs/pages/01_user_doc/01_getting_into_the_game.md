---
title: Getting into the game
---

## Getting into the game

### Getting the code

The NECI repository is stored on bitbucket. To gain access you need to
be invited. Contact one of the repository administrators [Simon Smart
([simondsmart@gmail.com](simondsmart@gmail.com)), George Booth
([george.booth24@gmail.com](george.booth24@gmail.com)) and Nick Blunt
([nsb37@cam.ac.uk](nsb37@cam.ac.uk))] who will invite you. If you
already have a bitbucket account let the repository administrators know
the email address associated with your account.

You will receive an invitation email. Please accept this invitation, and
create a bitbucket account as prompted if necessary.

To gain access to the NECI repository, an ssh key is required. This can
be generated on any linux machine using the command\footnote{asdfasdf}
[@dijkstra1959note]

```bash
ssh-keygen -t rsa -b 2048
```

This will create a private (`~/.ssh/id_rsa`) and a public key file
(`~/.ssh/id_rsa.pub`).

The private key must be kept private. On the bitbucket homepage, go to
account settings (accessible from the top-right of the main page), and
navigate to “SSH keys”. Click “Add key” and add the contents of the
public key. This will give you access to the repository.

You can now clone the code into a new directory using the command

```bash
git clone git@bitbucket.org:neci_developers/neci.git [target_dir]
```

### Required libraries

NECI requires some external software and library support to operate:

-   **MPI**<br>
    For builds of NECI intended to be run in parallel, an implementation
    of MPI is required. NECI has been heavily tested with OpenMPI, and
    MPICH2 and its derivatives (IBM MPI, Cray MPI, and Intel MPI).

-   **Linear algebra**<br>
    NECI makes use of the linear algebra routines normally contained in
    BLAS/LAPACK. There are a number of different packages which provide
    these routines, and are optimised for different compilers and
    platforms. NECI has been built and tested with either the AMD Core
    Math Library (ACML), the Intel Math Kernel Library (MKL), or the
    more general Basic Linear Algebra Subprograms (BLAS)/Linear Algebra
    Package (LAPACK) combination.

-   **HDF5 (optional)**<br>
    To make use of the structured HDF5 format for reading/writing
    POPSFILES (files storing the population of walkers, and other
    information, to restart calculations). This library should be built
    with MPI and fortran support (`--enable-parallel`
    `--enable-fortran` `--enable-fortran2003`).

-   **FFTW (optional)**<br>
    A small number of options in NECI, which are not enabled by default,
    require Fast Fourier Transforms. If these are re-enabled then the
    Fastest Fourier Transform in the West (FFTW3) library is required.

If combinations of these choices are made other than those most commonly
used then either the configuration files, or the resultant Makefile,
will need to be modified.

On the majority of machines available to the Alavi group and department,
the compilation environment is managed using the `module` command.
Documentation for that command is available from the command
`module help`. The most commonly used command is to load a module, using
the command

```bash
module load <module_name>
```

Installing and configuring the module system on private machines is far
beyond the scope of this document. Configuring your user account to use
modules may required modifications to your `.bashrc` file, depending on
the local machine configuration. Please contact the local IT
administrators or Simon Smart for further advice.

A number of standard combinations of modules present themselves. Where
an asterisk is presented, any version of the module can be used. Where
the version is under specified, the latest module should be used. This
will occur by default.

-   **gfortran (Cambridge)**<br>
    `mpi/mpich2/gnu/1.4.1p1 acml/64/gfortran/*/up`

-   **ifort (Cambridge)**<br>
    `ifort/64 mpi/openmpi/64/intel12 acml/64/ifort/*/up`

-   **PGI (Cambridge)**<br>
    `pgi/64 mpi/openmpi/64/pgi12 lapack/64/pgi` The lapack module might
    not be available to all users. Please contact Simon Smart if
    required.

-   **ifort (Max Planck FKF)**<br>
    `ifort mpi.intel` Note that the MKL library is included in the ifort
    module.

-   **ifort (hydra)**<br>
    `git intel mkl mpi.ibm hdf5-mpi cmake`

### Building the code (using CMake)

There are two ways of building NECI. The recommended approach is to use
the cmake build system. For legacy purposes, and for more explicitly
customised build configurations, the older Makefile system may also be
used.

CMake allows building the code in a separate directory. One directory
should be used per configuration that is to be built. This module can be
a subdirectory of the NECI directory, or otherwise. With the `build`
directory as the current working directory, execute

```bash
cmake [-DCMAKE_BUILD_TYPE=<type>] <path_to_neci>
```

pointing CMake at the root directory of the cloned NECI repository.

At this point CMake will automatically configure NECI according to the
currently loaded modules, and available compilers and libraries. If a
different set of modules or compilers are to be used a fresh directory
should be initialised (or equivalently all contents of the directory
deleted).

By default CMake will configure a Release build. If this is not desired
an alternative build type may be specified by setting the
`CMAKE_BUILD_TYPE` above. The available options are `Debug` (no
optimisations, all checking enabled), `Release` (optimisations enabled),
`RelWithDebInfo` (same as release, but with debugging symbols retained)
or `Cluster` (inter-procedural optimisations enabled, with very long
compile times).

@note
If there is an available HDF5 library, which is compiled with support
for MPI and for the Fortran compiler in use, then CMake will happily
make use of it. Otherwise support for HDF5 POPSFILES will be disabled by
default.

The CMake configuration for NECI contains the functionality to download,
compile and use HDF5. To do this run CMake with the
`-DENABLE_BUILD_HDF5=ON`.<br>
@endnote

The code is then built using the command

```bash
make [-j [n]] [neci|kneci|dneci|mneci]
```

The optional flag `-j` specifies that the build should be performed in
parallel (with up to an optional number, `n`, threads).

The final argument specifies whether the normal code (`neci`) should be
built, the complex code (`kneci`) the double run code (`dneci`) or the
multi run code (`mneci`) are built. If not specified, all targets are
built.

#### Options

Whilst every effort has been made to provide NECI with sensible default
options, the user may wish to play around further. To (dis)able an
option, the following should be passed as an argument to cmake:

```bash
-DENABLE_<option>=<(ON|OFF)>
```

The following options are available. Where an option is default "on", if
the required libraries are not available, the option will be disabled
and this will be noted in the build summary.

-   **BUILD\_HDF5**<br>
    Build the hdf5 library from source, and use that instead of one
    provided by the system.

-   **HDF5**<br>
    Make use of hdf5 for popsfiles (default=on).

-   **FFTW**<br>
    Functionality requiring FFTW (default=on).

-   **MOLCAS**<br>
    Build with the \_MOLCAS\_ flag (default=off).

-   **MPI**<br>
    Build with parallel functionality (default=off).

-   **SHARED\_MEMORY**<br>
    Use shared memory for storing integrals (default=on).

-   **WARNINGS**<br>
    Compile with verbose compiler warnings (default=off).

### Overriding configuration options

One of the aims of the CMake tool is to make build configuration as
black-box as possible. The build system should normally detect the
compilers in use automatically. The detected compilers may not, however,
be the ones desired, or the build system may fail to find functionality
that exists. The complier to use can be overridden by arguments passed
to CMake:

```bash
cmake -DCMAKE_<lang>_COMPILER=XXX <neci_dir>
```

where `<lang>` may be `Fortran`, `CXX` or `C` as appropriate. This
should specify the command to use which may be a compiler available in
the environmental `PATH`, or and absolute path to the compiler to use.

Once CMake has determined the compiler to use, it determines the
compilation and linker flags automatically. A number of overrides have
been definied for NECI. These may be found in the `cmake/compiler_flags`
director, and are segregated by files named according to both the
compiler vendor and the language.

Compiler flags are added in a specific order. The flags defined in the
`NECI_<lang>_FLAGS` variable are applied to all builds, whereas those in
`NECI_<lang>_FLAGS_<type>` are only applied to builds with the
appropriate build type (with the exception of when type is equal to
`CLUSTER` when the flags are appended to those used in `RELEASE` mode to
enable extra inter-file optimisations).

There are also flags to control how things are linked, what options are
passed to the compiler to enable compiler warnings, and flags that
depend on 32 or 64 bit builds. These should be self explanatory from
reading the files in `cmake/compiler_flags`.

If these defaults are insufficient, the compilation flags may also be
overridden. Any arguments passed to cmake of the form

```bash
cmake -DFORCE_<lang>_FLAGS[_<type>]
```

override the corresponding `NECI_<lang>_FLAGS[_<type>]` flags.
Essentially any `NECI_*` flag may be overridden on the command line with
a `FORCE_*` flag (although it is possible that some of these have been
accidentally ommitted from the implmentation).

#### Toolchain files

Overriding all of the CMake variables on the command line is cumbersome
and error prone. Various sets of overrides can be combined into a
toolchain file, which can be passed to CMake:

```bash
cmake -DCMAKE_TOOLCHAIN_FILE=<toolchain_file> <neci_dir>
```

These toolchain files can specify the entire chain of compilers, flags
and libraries if desired. For examples see the toolchains/ directory in
the NECI repository.

Additionally to the flags described above, the `CMakeForceCompiler`
functionality may be used (over and above just setting the
`CMAKE_<lang>_COMPILER` variable). These macros entirely disable the
autodetection of compiler properties within CMake (per language), and
will require all flags that are not in the `cmake/compiler_flags`
directories to be specified manually. As an example:

```cmake
include(CMakeForceCompiler)

CMAKE_FORCE_C_COMPILER       ( gcc GNU )
CMAKE_FORCE_CXX_COMPILER     ( g++ GNU )
CMAKE_FORCE_Fortran_COMPILER ( mpif90 GNU )
```

This will force the use of the commands `gcc`, `g++`, and `mpif90`, and
will set the `CMAKE_<lang>_COMPILER_ID` variable to `GNU` such that the
compiler flags set in the `cmake/compiler_flags/GNU_*.cmake` files are
used. Because this turns off any autodetection, CMake will not
autodetect that the c++ standard library needs to be linked in to
combine c++ and Fortran files. This will need to be corrected manually,
using:

```cmake
set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ )
```

Further internally required libraries may be required. In a normal
build, these are output in the build summary under "Implicit C++ linker
flags". A comprehensive documented example is found in
`toolchains/gfortran-openmpi.cmake`.

**Archer**

As an example, we can consider the Archer supercomputer, located at EPCC
in Edinburgh. This machine uses compiler wrapper scripts (written by
Cray) to provide much of the functionality automatically, but this
defeats the CMake auto-configuration system. To build on archer the
cmake command

```bash
cmake -DCMAKE_TOOLCHAIN_FILE=<neci_dir>/toolchains/archer.cmake <neci_dir>
```

(As we already know about Archer, we autodetect that you are running on
it, and CMake will fail with a message containing these instructions).
used, then the build will fail.

#### Compilation on ADA

Begin by clearing out the build environment. At the terminal execute:

```bash
module purge
module load environments/addons/cmake-2.8.2
```

Next, for GNU:

```bash
module load environments/programming/gcc-4.8.2
```

Or Intel:

```bash
module load compilers/intel/15.0.0.090
module load mpi/openmpi/1.8.2/intel15.0-threads
```

then run Cmake as normal.

#### Overriding packages required for options

There are a number of packages that are required to build NECI (such as
something providing a LAPACK-like interface), or which are required to
enable certain options (librt is required to enable shared memory).

If the package searching fails, there are a number of variables that can
be set on the CMake command line, or in a toolchain file as appropriate:

-   **NECI\_FIND\_&lt;package&gt;**<br>
    If this is set to OFF, the package searcher will not be executed,
    and the associated option will be enabled without disabling the
    option.

-   **&lt;package&gt;\_FOUND**<br>
    For many of the package searchers, this disables the searching from
    inside the package rather than outside. In general, the first option
    is preferred.

-   **&lt;package&gt;\_LIBRARIES**<br>
    The libraries to be linked in for use of the specifide package. If
    package searching is disabled, then to use the package this needs to
    be filled in explicitly.

-   **&lt;package&gt;\_DEFINITIONS**<br>
    Any additional compiler flags that are required to use the package.

-   **&lt;package&gt;\_INCLUDE\_PATH**<br>
    The location of any C or C++ header files, or fortran module files
    to be used during compilation

NECI has a some special package finders, called `MPI_NECI`,
`LAPACK_NECI` and `HDF5_NECI`. To a large extent they are just wrappers
around the underlying finders provided with CMake, but they implement
some additional logic, such as automatically substituting MKL for
LAPACK, checking the MPI compiler being used, and providing the capacity
to build hdf5 in the source tree.

As a result, when overriding these packages, `MPI_NECI`, `LAPACK_NECI`
and `HDF5_NECI` should be substituted for `<package>` above.

### Configuring builds (Makefile system)

A specific configuration for building NECI is initialised by using the
command

```bash
./tools/mkconfig.py config_name [-g]
```

The configuration names correspond to the configuration files contained
in the config directory. If the flag `-g` is used, then a debug
configuration will be created, and otherwise an optimised one.

There are a number of different configurations for differing systems,
library and compiler setups. The following are some of the more
important, and most likely to be used.

-   **gfortran\_simple, ifort\_simple**<br>
    These are the basic configurations, set up for the gfortran and
    ifort compilers. For development these are the most likely
    configurations to be used. On personal machines the easiest
    environment to install is gfortran and openMPI, using the
    `gfortran\_simple` - see note regarding libraries below.

-   **fkf\_ifort**<br>
    The MPI and ifort installation at the FKF in Stuttgart is different
    to that available in most locations. Use this config file to compile
    there.

-   **PC-ifort64-MPI-TARDIS, PC-ifort64-MPI-HYDRA**<br>
    The ifort compiler supports additional (expensive) optimisations
    during the link stage. For production runs on supercomputers these
    should be used. The -TARDIS config file is the normal one to use,
    with the -HYDRA for the differing environment available on HYDRA.

To specify a default configuration file to use on a particular machine,
create a symbolic link called `.default` in the config directory to the
appropriate configuration file. This allows `mkconfig.py` to be
usedwithout specifying the configuration name.

To compile NECI, a number of linear algebra routines are required. The
relevant routines are available in the AMD Core Math Library (ACML), the
Intel Math Kernel Library (MKL) and in the more widely available
combination of Basic Linear Algebra Subprograms (BLAS) and the Linear
Algebra Package (LAPACK). The configuration files above make some
assumptions about which packages are available, which are not always
correct.

In particular, the elements of the linker lines

```bash
-lacml
-lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
-lblas -llapack
```

are in principle interchangeable. On personal development machines it is
easiest to install BLAS and LAPACK, but these are generally less
performant so are not used by default. These lines can be substituted in
the generated `Makefile` before compilation.

### Building the code (Makefile system)

The code is built using the command

```bash
make [-j [n]] [neci.x|kneci.x|dneci.x|mneci.x|both|all]
```

The optional flag `-j` specifies that the build should be performed in
parallel (with up to an optional number, `n`, threads).

The final argument specifies whether the normal code (`neci.x`) should
be built, the complex code (`kneci.x`), the double run code (`dneci.x`),
the multiple run code (`mneci.x`) or both the normal and complex codes
(`both`) or all of the above and various extra utilities (`all`).

### Git overview

It is essential if you plan to do developmental work to get familiar
with the source-code management software ‘git’. The code will get
unusable exponentially quickly if all development and new ideas are
hacked into the master branch of the code. The nature of research is
that most things probably won’t work, but you want to implement them and
test relatively quickly, without requiring a standard of code that will
remain usable in perpetuity. To avoid an inexorable increase in code
‘clutter’, it is essential to work in ‘branches’ off the main code. For
a more detailed introduction to the git package, see
[git-scm.com/book/en/v2/getting-started-git-basics](git-scm.com/book/en/v2/getting-started-git-basics).
In short, the workflow should be:

1.  Branch off a clean master version to implement something

2.  Test and develop in the branch

3.  Regularly merge the new code from the master branch into your
    personal development branch

4.  Once satisfied with the development, and that it is an improvement
    in scope or efficiency of the existing code, ensure it is tidy,
    commented, documented, as bug-free as possible, and tests added to
    the test suite for it. This may involve reimplementing it from a
    clean version of master if it can be done more efficiently

5.  Merge code back into master branch

A few potentially useful git commands in roughly the workflow described
above:

-   **git branch**<br>
    See what branch I am on. -a flag for all (inc. remote) branches.

-   **git pull origin master**<br>
    Update the master branch into the current local repository

-   **git checkout -b newbranchname**<br>
    Fork off current branch to a new branch called ‘newbranchname’

-   **git commit -a -m ‘Commit message’**<br>
    Commit a set of changes for the current branch to your local
    repository.

-   **git push origin branchname**<br>
    Push your current local branch called branchname to a new remote
    branch of the same name to allow access to others and secure storage
    of the work

-   **git checkout -b newbranchname –track origin/remotebranch**<br>
    Check out a branch stored on the remote repository, and allow
    pushing and pulling from the remote repository for that branch.

-   **git push**<br>
    Push the current branch to the remote branch that it is tracking.

-   **git merge master**<br>
    Merge the recent changes in master into your local branch (requires
    a pull first)

-   **git checkout master**<br>
    Switch branches to the master branch

-   **git merge newbranch**<br>
    Merge your code in ‘newbranch’ into your current branch (potentially
    master)

Each commit should contain one logical idea and the commit message
should clearly describe *everything* that is done in that commit. It is
fine for one commit to only contain a very minor change. Try and commit
regularly and avoid large commits. It is also a good idea to make sure
that code compiles before commiting. This helps catch errors that you
may be introducing and also allows the use of debugging tools such as
git bisect.

It should be noted that the ‘stable’ branch of the code, automatically
merged into from master upon successful completion of nightly tests, is
hosted on github on a public repository, and also pushed to the molpro
source code. The molpro developers will quickly send us angry emails if
poor code gets pushed into it from NECI, and I will be sure to forward
complaints onto the relevant parties!
