This directory contains a collection of NECI test jobs. The test directories
contain benchmark files which are the output of 'correct' NECI calculations on
the corresponding tests. These benchmarks are used to check for errors
introduced into more recent versions of NECI.

The complete test suite (or any subset of the tests) can be run using testcode.
This is a python module which can run multiple tests and check the results
against the corresponding benchmarks.

testcode can be obtained from github:

$ git clone https://github.com/jsspencer/testcode ~/testcode2

Directory structure
===================

Important files in this directory are jobconfig, which contains details of
the tests (such as the number of processes to run with, which executables to use
and non-default tolerances in test comparisons) and userconfig, which contains
details of how to run jobs (such as which executables to use and how to specify
input files) and which benchmarks to compare against (specified by a benchmark
ID, or collection of benchmark IDs).

./bases contains integral files which are linked to by multiple tests.

./tools contains useful scripts, and extract.py, which is used by testcode to
extract the benchmarks values used to test against.

./neci, ./dneci and ./kneci contain tests which use the executables neci.x,
dneci.x and kneci.x respectively.

Running the test suite
======================

To run the complete set of tests, first compile neci.x, dneci.x and kneci.x,
i.e., in the directory above this one,

$ make neci.x && make dneci.x && make kneci.x

which are the three executables (currently) used by the test suite.

You can then run them all by running testcode in this directory:

$ ~/testcode2/bin/testcode.py

This will just tell you whether each test has passed, failed or been skipped.
To get more information you can increase the verbosity multiple times:

$ ~/testcode2/bin/testcode.py -v

or

$ ~/testcode2/bin/testcode.py -vv

which will give more information about failed tests.

You can run a subset of tests with the '-c' flag, for 'category'. For example,

$ ~/testcode2/bin/testcode.py -c neci

will run all tests in the neci directory (for which you only need neci.x
compiled, there is no need to also compile dneci.x and kneci.x in this case).

You can also run single tests, for example:

$ ~/testcode2/bin/testcode.py -c neci/serial/Ne_SS_Trial_Pops

Adding a new test
=================

To add a new test, first create a directory in the appropriate subdirectory
(for a parallel neci.x job, inside ./neci/parallel). In this directory add the
test's input file (with a name ending in '.inp') and any necessary additional
files, such as integral files or POPSFILEs.

You should then add these files to git, for example:

$ git add neci/parallel/new_test

If the test is added to a new subdirectory then you may need to add it to
./jobconfig. If you added it to an already-exisiting directory, such as
neci/parallel, then it should be automatically found by testcode using the
globbing in jobconfig.

You must then create a benchmark file by running the test suite. To create a
new benchmark for *only* the new test, run

$ testcode2/bin/testcode.py make-benchmarks -ic neci/parallel/new_test

This should run just the new test. You will be told that the test has failed,
and asked if you would like to set the new benchmark. If you believe that the
test has run correctly then do so.

'-i' tells testcode to 'insert' the new benchmark ID at the start of the old
list of benchmarks (located in ./userconfig). When testcode is run later, it
will use the benchmark files with these IDs to compare against.

If you don't include '-i' then testcode will remove all previous benchmarks.
This is useful if you want to reset benchmarks for the whole test suite, which
can be done with:

$ testcode2/bin/testcode.py make-benchmarks

Finally, once this is done you need to add the new benchmark file to git.
If the new benchmark ID is, for example,

dneci-c3462e0.kneci-c3462e0.neci-c3462e0

then you can do:

$ git add neci/parallel/new_test/*dneci-c3462e0.kneci-c3462e0.neci-c3462e0*

You should then commit the newly added files. You might like to try running
testcode on the new test, and making sure it runs and passes as expected,
confirming the that the test was added correctly:

$ ~/testcode2/bin/testcode.py -c neci/parallel/new_test
