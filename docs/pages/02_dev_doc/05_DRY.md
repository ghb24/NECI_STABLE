---
title: Don’t Repeat Yourself (DRY)
---

## Don’t Repeat Yourself (DRY)

When information becomes duplication in software, eventually the people
that knew about the duplication will forget. And then the information
will be changed — but at least one of the duplicates won’t be. This
introduces bugs that are extremely difficult to track down.

The Don’t Repeat Yourself (DRY) principle is one that say a developer
should systematically, and always, avoid duplication of information.

NECI is a terrible example of this, but it has been improved over time
with a lot of effort.

Information is a very broad term. There are many types of duplication
that can occur. A non exhaustive list of some types of duplication (and
what can be done about them) follows.

-   **Algorithm duplication across data types**<br>
    There are many algorithms that are either the same, or similar,
    across many different data types. The logic involved in these should
    be written once.

    Major examples are the sorting routines in `sort_mod`, general
    utilities in `util_mod`, shared memory in `shared_alloc` and MPI
    routines in `Parallel_neci`. Prior to implementing a generalised
    quicksort there were 37 different sort routines in NECI, using
    different sort methods, and containing different bugs.

    This duplication should be controlled using templating, as described
    in
    section <a href="#sect:templating" data-reference-type="ref" data-reference="sect:templating">1.7</a>.

-   **Logic duplication across source files**<br>
    If the same chain of decision making is recurring in different
    regions of the code, these should be abstracted into their own
    subroutine which is called from each location. This prevents the
    logic being duplicated, and then diverging.

    Numerous bad examples of this still persist in NECI.

-   **Duplication of data**<br>
    Compile time constants should only be specified in one place. A
    large proportion of these are found in the `lib/cons_neci.F90`
    source file. Other examples include the layout of the bit
    representations (`BitReps.F90`). A significant proportion of these
    vary depending on the compile configuration, and prior to collecting
    them here the code was extremely fragile.

    Ongoing cases which are problematic include the use of the literal
    constant \(6\) to specify output to stdout in statements such as
    `write(6,*)`, which doesn’t interact well with `molpro`.

-   **Duplication of representations in memory**<br>
    It is important to have a well defined canonical representation of
    data in memory. The same data should not be allowed to become
    duplicated in multiple places.

    Temporary arrays, with working data copied into them, should be
    clearly temporary and discarded as soon as not necessary. If the
    primary data shifts to a new location, the old storage should be
    deallocated (if possible), or damaged so that attempts to use it
    fail loudly (such as putting a value of \(-1\) into a variable that
    would normally hold an index).

    Avoid situations where code might work by accident.

    An ongoing situation of this type is the array variable `nBasisMax`.
    It shadowed a large number of global control variables, and there
    are still locations in the code where its value is used in
    preference to the global control value as these values diverge, and
    its value is the one that works in some obsolete code.

There is one, major, exception to this rule. Code that only exists for
testing purposes (such as the contents of `ASSERT` statements, or unit
tests) may be as explicitly duplicated as desired. Their *purpose* is to
explicitly flag up when anything elsewhere changes - so the risk of them
getting out of sync with the code base is their purpose.

### Procedure pointers (function pointers)

One issue that becomes immediately obvious when repeated logic has been
abstracted into specific functions is that conditional logic is executed
*every* time certain actions are taken.

In many cases this is not very important, as it occurs high up the call
hierarchy, and the controlled code consumes the vast majority of the
execution time. However, the closer we get to the inner most tight
loops, the more expensive repeated conditional logic becomes. This is
particularly frustrating if the decision making is based on global
control parameters, and thus always results in the same code path being
taken in a simulation. If the decision lies against the branch
prediction metrics, then this is especially bad.

The canonical example of this is accessing the 4-index integrals, which
is performed very frequently.

In these case, it is a good idea to separate the decision making logic
from the execution, so that the conditional logic is only executed at
runtime. This can be done using *procedure pointers* (called function
pointers, or similarly functors in other programming languages.

These require defining the “shape” of a function call (i.e. its
arguments, and return values) in an `abstract interface`. A variable can
then be set to point at which of a range of functions with this
signature should be executed.

In NECI the global controlling procedure pointers are located in the
module `procedure_pointer`, which contains both the abstract interface
definitions, and the actual pointer variables. These variables can then
be used as functions throughout the code.

These procedure pointers are largely initialised in the routine
`init_fcimc_fn_pointers`, where decisions are made between types of
excitation generator, matrix element evaluation, etc. The procedure
pointers involved in integral evaluation are set in
`init_getumatel_fn_pointers`.

The use of procedure pointers generates a strict Fortran 2003 dependency
for NECI. We used to make use of a hacky abuse of the linker and
templating system to implement function pointers without language
support, but this was deprecated once compiler support for procedure
pointers was reasonably widespread.
