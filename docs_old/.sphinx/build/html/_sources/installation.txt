.. _installation:

============
Installation
============

Requirements:

* git.
       The NECI codebase, documentation and test suite are currently distributed
       only from git repositories hosted on CUC3 workstations (account required).
* subversion.
       Legacy code and a modified version of CPMD are distributed via the wwmm
       subversion repository (account required), which is kindly hosted by the
       Unilever Centre.
* LAPACK.
* BLAS.
* FFTW 3.x.

--------
Download
--------

NECI
----

The NECI source code can be checked out using::

.. code-block:: bash

   $ git clone scepter:/home/ajwt3/NECI.git

Various branches also exist, but they may or may not be stable or under
active development.

CPMD
----

Our modified version of CPMD (containing the necessary routines for
integration with NECI) can be downloaded from::

.. code-block:: bash

    $ svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/CPMD/branches/QMC/trunk CPMD

Again, there exist various branches and they may or may not be stable.

The default setting used in the compilation scripts assumes that the NECI
source is in a subdirectory of the CPMD source directory.  This can be changed
by using the command line options or setting them in the .compileconf file (see
below).

Note that the codebase is not currently under active development.

--------------
File structure
--------------

NECI
----

config/
    Directory containing configuration files used for building NECI on a
    variety of platforms and compilers.
docs/
    Documentation.
src/
    Source files.
tools/
    Tools which are useful for running and analysing NECI output.  These
    include a histogramming script and scripts which automate sending a
    **SOFTEXIT** signal to a running instance of NECI after a set amount
    of time.
utils/
    Useful scripts for compiling and working on the code base.

The following directories are created during compilation

dest/
    Compiled objects.
bin/
    Executables.
lib/
    NECI libraries for linking to VASP/NECI.

The major compiled files are:

NECI/bin/neci.x
    Standalone neci executable.
NECI/lib/neci-cpmd.a
     NECI library for integration with the CPMD gamma-point code.
NECI/lib/neci-vasp.a
     NECI library for integration with the CPMD k-point code.
NECI/lib/neci-cpmd.a
    NECI library for integration with the VASP gamma-point code.
NECI/lib/neci-vasp.a
    NECI library for integration with the VASP k-point code.

All of these are actually symbolic links.  The binaries and libraries are
unique to the optimisation level and configuration used.  The above files
merely point to the most recent appropriate binary or library. 

CPMD
----

CONFIGURE/
    Directory containing the configuration files for a variety of platforms.
dest/
    Build directory for gamma-point code.
kdest/
    Build directory for k-point code.
gcpmd.x
    Gamma-point executable of the CPMD-NECI code (links to CPMD/dest/cpmd.x).
    Must not be used for k-point calculations!
kcpmd.x
    k-point executable of the CPMD-NECI code (links to CPMD/kdest/cpmd.x).
    Must not be used for gamma-point calculations!

-----------
Compilation
-----------

The settings provided are mainly for CUC3 machines and other computer clusters
that the Alavi group has access to and these might need to be adjusted in order
to compile in different environments.  Most of the time only the paths and
flags for the FFTW, LAPACK and BLAS libraries, given in the LFLAGS variable,
will need to be adjusted (if that).  It is necessary to use the same compiler
to compile CPMD as was used to produce the NECI library for CPMD and similarly
for integration with VASP.

We have compiled and tested the codebases with the gfortran (4.2 and later),
Portland, Pathscale and Intel compilers in 32 and 64 bit.  

NECI
----

To compile NECI run::

.. code-block:: bash

    $ ./tools/mkconfig.py platform
    $ make

The platform is a filename in the config directory.  The mkconfig.py script
has some useful options.  Run::

.. code-block:: bash

    $ ./tools/mkconfig.py --help

to see information on them.

The objects are compiled to dest/platform/optimised/real (or complex, for
compiling the complex version for libraries).  The resultant executable,
neci.platform.optimised.x, is placed in the exe directory.

If the debug flag [-g] is given to mkconfig.py, then the debug configuration is used and 
the filenames and paths contain debug rather than optimised.

For convenience, bin/neci.x is a symbolic link to the most recently compiled
executable.

Compiling libraries works in much the same way.

There are several goals defined in the makefile.  Run

.. code-block:: bash

    $ make help

to see the most useful.

Configuration files use a simple ini format and adding a new configuration is
easy. See the supplied configurations for examples and/or comments at the start
of the mkconfig.py script.

CPMD
----

For historical reasons, CPMD and NECI are much more closely interwoven
than NECI is with VASP.

The repository version of CPMD depends upon NECI, thus NECI must be
compiled first.  Various scripts take care of this for us.

CPMD and NECI have CONFIGURE subdirectories.  Each file in a CONFIGURE
subdirector contains the necessary information to compile the code using
a certain compiler.

Please note that not all the CPMD configure scripts will work: many of
them were supplied with CPMD and have never been used with our modified
version.  Please only use platforms that exist in the NECI/config directory as
well as the CPMD/CONFIGURE directory.  It is easy to make your own
configuration files using existing ones as a template.

Note that some of the platforms obtain LAPACK and BLAS as part of atlas,
ACML or MKL, so this will need to be changed if different source libraries
are used.

Two versions of CPMD (and the corresponding NECI library) exist, gcpmd.x
(for gamma-point calculations) and kcpmd.x (for k-point calculations),
to take advantage of some substantial memory savings when running
gamma-point calculations, as all wavefunctions are then real.  This is
controlled via a C pre-processing statement.  As neci.x is only for
molecular calculations, it only exists in one form.

runmake.sh is a script for CPMD which controls the process of creating
makefiles and compiling the NECI/CPMD hybrid.  A "platform" is just the name of
one of the configuration files in the CONFIGURE subdirectories.

Run::

.. code-block:: bash

     $ ./runmake.sh -h

to see the various options and default behaviour of the runmake.sh script.

Note that runmake.sh produces new makefiles for CPMD **and** for NECI,
and compiles neci.x, the NECI libraries needed for CPMD, and gcpmd.x
and kcpmd.x.  To aid the speed of recompilation, the real and complex
objects are compiled into separate directories.

runmake.sh defaults to compiling the codebases using the Portland 64-bit
compiler, if a platform is not specified either via the command line or given
in .compileconf, a text file which contains the name of the desired platform.
Note that runmake.sh will use the same platform for both the CPMD and NECI
makefiles.

The CPMD source directory also contains a controlling Makefile
to further help the make process (and generally just acts as a wrapper
for the runmake.sh script).  Run::

 make help

in each directory to see the various targets available.

To quickest way to compile both CPMD and NECI is to run::

 [CPMD]$ make all

from within the CPMD source directory.

.compileconf
^^^^^^^^^^^^
The .compileconf file is not under source code management and allow local defaults
to be set.

When running runmake.sh, please note that it uses the CPMD .compileconf
information for compiling NECI, rather than the user's default NECI platform.
This is to ensure that the same platform is used for both.

The settings in .compileconf are overridden by command line options but override
any defaults in the runmake.sh script.

.compileconf in its simplest (and oldest) form simply contains the name of the desired
platform, e.g.::

    PC-ifort64

will use the PC-ifort64 platform as the default.

.compileconf can also be used to set local defaults for more variables---see
the comments in runmake.sh for more details.  Defaults are used for any
variables not set in the .compileconf file.

For example, to set different defaults for the platform and the location of the
NECI source, .compileconf in the CPMD directory would look like::

    platform=PC-ifort64
    NECIsrc=~/NECI/source

---------------
Developing NECI
---------------

Default platform
----------------

If no platform is given on the command line to the mkconfig.py script, then it
looks for a .default file in the config subdirectory, which is used if it exists.

It is thus simple to set a local default platform by creating a symbolic link from 
the desired platform file to config/default.

Working in the src subdirectory
-------------------------------

Note that any goal defined in the Makefile in the root NECI directory can also
be run from the src directory (see src/Makefile for how this is achieved).

Adding new files
----------------

Just create the new source file in the src subdirectory.  The makefile will
pick it up automatically.  Note that *all* .F, .F90, .c and .C files in src
will be compiled and linked to form neci.

Updating dependencies
---------------------

The dependency list is created automatically if it does not exist using.  If this
needs to be updated (e.g. due to a new source file or because development results in
additional dependencies), run::

.. code-block:: bash

    $ make depend

Note that this causes an infinite loop if make 3.80 or earlier is used.  As an alternative::

.. code-block:: bash

    $ make rmdeps

deletes all dependency files, which are then automatically regenerated the next time make
is run to produce any target.  This should be used on older systems.

Note that all platforms share the same set of dependency files.

Testing on different platforms
------------------------------

Because each neci executable is given a file name dependent upon the
configuration and optimisation level used, it is easy to test multiple
configurations (e.g. parallel and serial) within one local repository.

.. code-block:: bash

    $ ./tools/mkconfig.py -f Makefile.serial PC-ifort64
    $ ./tools/mkconfig.py -f Makefile.mpi PC-ifort64-MPI
    $ make -f Makefile.serial
    $ make -f Makefile.mpi

This will compile 2 executables.  Recompiling will only involve recompiling
any changed files. 

tags
----

Thanks to Alex for pointing this out.  It is useful to be able to work on the
source code both from the main directory and the src directory and have the
same tags file work in both cases.

You can create a tags file by running:

.. code-block:: bash

    $ make tags

and tell vim to search the current directory and parent directory (which
amounts to working in the main and src directories respectively for our case)
by doing::

    :set tags=./tags,../tags

within vim or placing::

    set tags=./tags,../tags

in your $HOME/.vimrc 

--------
testcode
--------

testcode is a set of scripts written by James Spencer that is used to
check that our programs produce the same results as they did before.
It is useful both for development work, to ensure that regression issues
are avoided, and testing successful compilations.

Every night the latest version of the codebase is checked out of the
subversion repository and tested against a variety of compilers, giving
confidence in the continued stability of the codebase.

The testcode scripts, NECI test suite and NECI/CPMD test suite can
be obtained respectively via::

.. code-block:: bash

    $ git clone scepter:/home/jss43/alavi_group/testcode.git
    $ git clone scepter:/home/jss43/alavi_group/testsuite_neci.git
    $ git clone scepter:/home/jss43/alavi_group/testsuite_cpmd.git

Please see the testcode documentation for more details.
