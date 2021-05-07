---
title: Developer's guide to NECI
author: - Simon Smart, Nick Blunt, Oskar Weser and George Booth
---

[TOC]
# Working in NECI

In many ways, NECI is a fairly hostile environment to code, especially
for inexperienced software developers, or those who are not familiar
with the ideosynchrasies of different versions of FORTRAN/Fortran.

In the following sections we aim to give some general guidance for
working in the NECI codebase, and

@note
As a general guideline, programs should be written to fail as loudly and
as early as possible. This pushes the job of finding errors and
debugging up the tree.

1.  Just “looks wrong” to the programmer.

2.  Syntax highlighting in the editor makes a mistake stand out.

3.  Error detected by compiler.

4.  Error detected by linker.

5.  Runtime error detected by debug sanity checks through code.

6.  Runtime error caused in code at a different location to the bug.

7.  Runtime error only occurs when running in parallel, or in bigger
    calculations.

8.  Simulation appears to run correctly, but gives obviously wrong
    output.

9.  Simulation appears to run correctly, but gives subtly wrong output.

We strongly aim to be at the top of the list.
@endnote

## Code conventions


If a function is `pure` and operates on scalars, add the `elemental`
attribute to automatically map it elementwise onto arrays.

The following toy function is pure, can be applied onto arrays and
scalars alike. If the exponent is ommited, it defaults to squaring.

    elemental function pow(x, n) result(res)
        integer, intent(in) :: x
        integer, intent(in), optional :: n
        integer :: n_

        integer :: i

        def_default(n_, n, 2)

        res = 1
        do i = 1, n_
            res = res * x
        end do
    end function

    pow([1, 3, 5]) -> [1, 9, 25]
    pow([1, 3, 5], 3) -> [1, 27, 125]

```Fortran
elemental function pow(x, n) result(res)
    integer, intent(in) :: x
    integer, intent(in), optional :: n
    integer :: n_

    integer :: i

    def_default(n_, n, 2)

    res = 1
    do i = 1, n_
        res = res * x
    end do
end function

pow([1, 3, 5]) -> [1, 9, 25]
pow([1, 3, 5], 3) -> [1, 27, 125]
```

**Data types**
With the exception of small integers being directly assigned to known
integer variables, or used in loop counters, all constants should have
their types explicitly specified. The available types are described in
section
<a href="#sect:cons-types" data-reference-type="ref" data-reference="sect:cons-types">2.1</a>.

Most specifically, the `*D*` specifier and the floating point type
`double precision` should never be used. Examples such as `1.D0` should
be replaced with `1.0_dp`, and the data types `real(dp)` and
`complex(dp)` should be used.

As compilers have moved between 16-bit, 32-bit and 64-bit, there is
ambiguity about whether `double precision` should mean a 32-bit, 64-bit
or 128-bit floating point value, depending on the age of the compiler
and which compiler is used. This can cause chaos and difficult to track
runtime bugs that appear only on certain machines.

For Hamiltonian matrix elements (that may be real or complex depending
on build configuration) the custom (preprocessor defined) data type
`HElement_t` should be used, which resolves to either `real(dp)` or
`complex(dp)`.



### Variable substitution

The templater has extremely powerful mechanisms to manipulate the types
of variables. Variable manipulation is enabled by using a key in the
key-value pair section that begins with `type`. In particular, the code
is able to manipulate the number of dimensions that different arrays
have.

As an example, take a routine which is passed an array, and a value that
could be an element of that array (such as is necessary for a binary
search), such that


```Fortran
subroutine example(arr, elem)
    %(type1)s :: arr(:)
    %(type1)s :: elem()
    ...
end subroutine
```

If `type1` is a scalar value, this does a substitution exactly as would
be expected:

<!-- TODO should have ⇒ where TODO is -->

```Fortran
[int]
type1 = integer(int32)
```
becomes

```Fortran
subroutine example(arr, elem)
    integer(int32) :: arr(:)
    integer(int32) :: elem()
    ...
end subroutine
```
However, it may be necessary for the value which is being considered in
the array to itself be an array. An example of this would be the bit
representations used in NECI — a list of of these is a two dimensional
array, and any intermediate values would be arrays themselves.

