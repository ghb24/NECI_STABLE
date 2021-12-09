---
title: Testing
---

## Testing

NECI comes with a set of tests in the test\_suite directory. Each of
these tests have a benchmark file. When you run the tests, the results
of your calculation will be compared against those from the benchmark
files. If the values of certain results agree to within a predefined
tolerance, the test will pass.

The test suite is run using a python program called testcode2. This
program will call the desired tests, compare the results against
benchmarks, and let the user know the outcome of each test.

You can clone testcode2 from github with the following command:

```bash
$ git clone https://github.com/jsspencer/testcode ~/testcode2}
```

There are tests for each of the three executables: neci, mneci and
kneci. These are stored in the three directories with the corresponding
names. These directories are then divided into further directories for
the different types of tests. For example, mneci has subdirectories
called rdm, excited\_state, kpfciqmc and so on, with tests for each of
these corresponding features of NECI.

To run the entire test suite, just do

```bash
$ ~/testcode2/bin/testcode.py
```

in the test\_suite directory (assuming you cloned testcode2 to your home
directory). To run this, you will have to have all of neci, mneci and
kneci compiled. However, you can also run a subset of tests. For
example, to run all mneci tests do

```bash
$ ~/testcode2/bin/testcode.py -c mneci
```

or to run a particular single test do

```bash
$ ~/testcode2/bin/testcode.py -c mneci/rdm/HeHe_int
```

or to run two particular tests do

```bash
$ ~/testcode2/bin/testcode.py -c mneci/rdm/HeHe_int -c mneci/rdm/HeHe_real
```

By default, testcode will just tell you whether or not the test passed.
If the test failed, you can get further information by increasing the
verbosity of the output. For example,

```bash
$ ~/testcode2/bin/testcode.py -c mneci/rdm/HeHe_int -v
```

or

```bash
$ ~/testcode2/bin/testcode.py -c mneci/rdm/HeHe_int -vv
```

which will tell you why the individual test did not pass. You can also
use the verbosity flags when running the entire set of all tests.

### Adding a new test

To add a new test, first create a directory in the appropriate
subdirectory (for a parallel neci job, inside ./neci/parallel). In this
directory add the test’s input file (with a name ending in ’.inp’) and
any necessary additional files, such as integral files or POPSFILEs.

You should then add these files to git, for example:

```bash
$ git add neci/parallel/new_test
```

If the test is added to a new subdirectory then you may need to add it
to ./jobconfig. If you added it to an already-exisiting directory, such
as neci/parallel, then it should be automatically found by testcode
using the globbing in jobconfig.

You must then create a benchmark file by running the test suite. To
create a new benchmark for *only* the new test then run, for example,

```bash
$ testcode2/bin/testcode.py make-benchmarks -ic neci/parallel/new_test
```

This should run just the new test. You will be told that the test has
failed, and asked if you would like to set the new benchmark. If you
believe that the test has run correctly then do so with `y`.

`-i` tells testcode to ’insert’ the new benchmark ID at the start of the
old list of benchmarks (located in ./userconfig). When testcode is run
later, it will use the benchmark files with these IDs to compare
against.

If you don’t include ’-i’ then testcode will remove all previous
benchmarks. This is useful if you want to reset benchmarks for the whole
test suite, which can be done with:

```bash
    $ testcode2/bin/testcode.py make-benchmarks
```

Finally, once this is done you need to add the new benchmark file to
git. If the new benchmark ID is, for example,

```bash
mneci-c3462e0.kneci-c3462e0.neci-c3462e0
```

then you can do:

```bash
$ git add neci/parallel/new_test/*dneci-c3462e0.kneci-c3462e0.neci-c3462e0*
```

You should then commit the newly added files. You might like to try
running testcode on the new test, and making sure it runs and passes as
expected, confirming the that the test was added correctly:

```bash
$ ~/testcode2/bin/testcode.py -c neci/parallel/new_test
```

### Unit tests

Unit tests should test discrete, small, elements of functionality.
Ideally these should be the smallest elements such that functionality is
then composed from “units” that have all been tested. By considering
each of the elements explicitly, it is possible to writ tests for
edge-case behaviour, where it is difficult to ensure that this behaviour
will be tested in a full integration test case.

Unit tests are found in the `unit_tests` directory. They are arranged in
subdirectories, each of which corresponds to one of the *files* of NECI
source code. There are a number of technical steps to integrating new
unit tests.

1.  A passing test is an executable that returns 0, and any other return
    value indicates failure. The developer has an entirely free choice
    to determine how they wish to write test executables.

2.  The library FRUIT is provided to assist in writing unit tests.
    Within a given test (i.e. testing a given function, or unit), there
    should be a *suite* of tests to cover all possbile cases. FRUIT
    provides helper functionality to keep track of which element of a
    suite is currently running (using the `TEST()` macros), check values
    (using `call assert_equals`, `call assert_true` and so forth), keep
    track of where errors occurred and report them in an easily readible
    form. The module can be imported with `use fruit`. See existing
    tests for examples.

3.  If the current directory of tests does not contain a
    `CMakeLists.txt` file, create it. Ensure that it contains a
    `foreach()` loop over the available
    `\${PROJECT_NAME}_CONFIGURATIONS` so that all the build
    configurations of neci get tested. The directory should be added to
    the main `CMakeLists.txt` in the `unit_tests` directory using the
    `add_subdirectory()` command.

4.  Add the test to the `CMakeLists.txt` file in the directory in which
    it resides using the `neci_add_test` command. This will require you
    to specify a name for the test, and the appropriate .F90 file.

5.  To test the appropriate configuration of neci, add `lib{k,m,d,}neci`
    to the LIBS line. To use the FRUIT helpers, add `fruit` to this
    line.

The easiest way to get these details correct is to copy existing
examples.

Unit tests can be run using the command

```bash
ctest [-R <regex>]
```

which is a built in part of the CMake toolkit. By default this will
execute all avaialable tests and print a report on successes and
failures. The optional `-R` flag specifies a regular expression, and
only tests matching this will be executed. There are a number of further
options available.

A test can also be run directly by running its executable manually. This
can be more straightforward for capturing the output of a failing test
whilst debugging. The executables are found in the same position in the
build directory as the source files are in the main repository
(i.e. `unit_tests/det_bit_ops/test_countbits.F90` gives
`<build_dir>/unit_tests/det_bit_ops/test_countbits`).
