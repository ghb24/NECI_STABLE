.. _installation:

============
Installation
============

Requirements:

* subversion
       The Alavi group currently only distributes code via the wwmm
       subversion repository (account required), which is kindly hosted
       by the Unilever Centre.
* LAPACK
* BLAS.
* FFTW 3.x.

--------
Download
--------

NECI
----

The NECI source code can be downloaded from::

 svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/NECI/trunk NECI

Various branches also exist, but they may or may not be stable or under
active development.

CPMD
----

Our modified version of CPMD (containing the necessary routines for
integration with NECI) can be downloaded from::

 svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/CPMD/branches/QMC/trunk CPMD

Again, there exist various branches and they may or may not be stable
or under active development.
The default setting used in the compilation scripts are to have the NECI source as
a subdirectory of the CPMD source.  This can be changed by using the
command line options or setting them in the .compileconf file (see below).
It is suggested that you place NECI as a subdirectory of the CPMD source,
as that is how the CPMD configuration files in the repository are set up.

VASP
----

Only users on the access list can download VASP from our repository::

 svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/VASP/vasp.4.lib vasp.4.lib
 svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/VASP/vasp.5 vasp.5

-----------
Compilation
-----------

CPMD and NECI
-------------

For historical reasons, CPMD and NECI are much more closely interwoven
than NECI is with VASP.  In addition, the tools used to compile CPMD and
NECI are very similar (but then they were written by the same people!).

The repository version of CPMD depends upon NECI, thus NECI must be
compiled first.  Various scripts take care of this for us.

CPMD and NECI have CONFIGURE subdirectories.  Each file in a CONFIGURE
subdirector contains the necessary information to compile the code using
a certain compiler.  Most of the time only the paths and flags for the
FFTW, LAPACK and BLAS libraries, given in the LFLAGS variable, will need
to be adjusted (if that).  It is necessary to use the same compiler to
compile CPMD and NECI.  We have compiled and tested the codebases with the
gfortran (4.2 and later), Portland and Intel compilers in 32 and 64 bit.
The mkconfig.sh scripts in the CPMD and NECI directories produce the
relevant makefiles, but this is better done via helper scripts.

Please note that not all the CPMD configure scripts will work: many of
them were supplied with CPMD and have never been used with our modified
version.  Please only use platforms that exist in the NECI/CONFIGURE
directory as well as the CPMD/CONFIGURE directory.  It is easy to make
your own configuration files using existing ones as a template.

Note that some of the platforms obtain LAPACK and BLAS as part of atlas,
ACML or MKL, so this will need to bechanged if different source libraries
are used.

Two versions of CPMD (and the corresponding NECI library) exist, gcpmd.x
(for gamma-point calculations) and kcpmd.x (for k-point calculations),
to take advantage of some substantial memory savings when running
gamma-point calculations, as all wavefunctions are then real.  This is
controlled via a C pre-processing statement.  As neci.x is only for
molecular calculations, it only exists in one form.

runmake.sh and compile are scripts for CPMD and NECI respectively which
control the process of creating makefiles and compiling the codebases.
They both refer to a "platform", which is just a name of one of the
configuration files in the CONFIGURE subdirectories.

compile by default generates new makefile (for both gamma-point and
k-point compilations) and does a clean build of neci.x::

 [NECI]$ ./compile -h
 usage: ./compile [options] [platform]
 Generate new makefiles and do a clean build of neci.x using platform as the configuration.

 If platform is not specified, then the platform given in .compileconf is used.
 If .compileconf also doesn't exist, then the default (PC-PGI64) is used.

 Options:
 -d Compile with the compiler debug options on.
 -f Fast: don't do a make clean before compiling.
 -m Only make new makefiles.

In contrast, runmake.sh has a different default behaviour, in that it
doesn't produce new makefiles by default, and does not do clean builds::

 [CPMD]$ ./runmake.sh -h
 Usage: ./runmake.sh [-c] [-h]
 Compile NECI (neci.x) and CPMD/NECI code for Gamma point (gcpmd.x)
 code and for k-point sampling (kcpmd.x).
 Warning: the option to set the NECI source directory are
 *only* used when a new makefile is produced (i.e. requires the -m or -p
 flag).

 Options:
     -c  Recompile only CPMD routines.
     -d  Generate new makefiles for debugging.  Recompile (at least)
         CPMD/qmc routines and all of NECI.
     -g  Compile only Gamma point code.
     -k  Compile only K-point code.
     -m  Generate new Makefiles by running mkconfig scripts in NECI and
         CPMD directories, using the default platform (PC-PGI64 unless 
         otherwise specified in .compileconf), and compiles.
     -n  Recompile only NECI routines.
     -p [32,64,platform] 
         Produce makefile for [pgi-32bit,pgi-64bit,platform]
         compilation, where platform is an alternative configuration (eg
         for gfortran).
     -s [NECI source directory]
         Set the location of the directory containing the NECI source code.
         Warning: must be used only when new makefiles are produced (i.e.
         when -m or -p are specified).
     -h  Print this message.

Note that runmake.sh produces new makefiles for CPMD **and** for NECI,
and compiles neci.x, the NECI libraries needed for CPMD, and gcpmd.x
and kcpmd.x.  To aid compilation, the dest subdirectories in the CPMD and
NECI source directories contains the compiled objects for the gamma-point
code and the kdest subdirectories contain the compiled objects for the
k-point code.

Both the NECI and CPMD scripts default to compiling the codebases using
the Portland 64-bit compiler, if a platform is not specified either via
the command line or given in .compileconf, which is a text file which
contains the name of the desired platform.  Note that runmake.sh will
use the same platform for both the CPMD and NECI makefiles.

The CPMD and NECI source directories also contain controlling Makefiles
to further help the make process (and generally just act as wrappers
for the runmake.sh and compile scripts).  Run::

 make help

in each directory to see the various targets available.

To quickest way to compile both CPMD and NECI is to run::

 [CPMD]$ make all

from within the CPMD source directory.

.compileconf
^^^^^^^^^^^^
The .compileconf files are not under source code management and allow local defaults
to be set.  It is used both in the CPMD and NECI compilation scripts.

When running runmake.sh, please note that it uses the CPMD .compileconf information
for compiling NECI, rather than the NECI .compileconf file.  This is to ensure that
the same platform is used for both.

The settings in .compileconf are overridden by command line options, but override
any defaults in the compilation scripts.

.compileconf in its simplest form simply contains the name of the desired
platform, e.g.::

    PC-ifort64

will use the PC-ifort64 platform as the default.

The CPMD .compileconf can be used to set local defaults for more
variables---see the comments in runmake.sh for more details.  For example, to set
different defaults for the platform and the location of the NECI source::

    platform=PC-ifort64
    NECIsrc=~/NECI/source

File structure
^^^^^^^^^^^^^^
NECI files:

**NECI/neci.x**
 Standalone neci-executable (links to NECI/dest/neci.x).
**NECI/dest/neci-cpmd.a**
 NECI library for CPMD gamma-point code.
**NECI/kdest/neci-vasp.a**
 NECI library for CPMD k-point code.
**NECI/dest/neci-cpmd.a**
 NECI library for VASP gamma-point code.
**NECI/kdest/neci-vasp.a**
 NECI library for VASP k-point code.

CPMD files:

**CPMD/gcpmd.x**
 Gamma-point executable of the CPMD-NECI code (links to CPMD/dest/cpmd.x).
 Must not be used for k-point calculations!
**CPMD/kcpmd.x**
 k-point executable of the CPMD-NECI code (links to CPMD/dest/cpmd.x).
 Must not be used for gamma-point calculations!

VASP
----

James has managed it.  It's not completely pleasant.  More to follow
once it's been made easier!

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

testcode and the set of test jobs (both for NECI and CPMD-NECI), can be
checked out of the subversion repository:

.. code-block:: bash

    svn checkout https://wwmm.ch.cam.ac.uk/svn2/groups/alavi/testcode testcode

Please see the testcode documentation for more details.
