---
title: Interfacing C (and C++) code
---

## Interfacing C (and C++) code

Developers have always mixed Fortran code with external routines
implemented in other languages, especially C. This has generally been
done on an ad-hoc basis by exploiting the naivety of the linker — in
particular that the linker will resolve dependencies with any symbol of
the specified name. This is useful, but introduces a number of potential
problems:

-   **Name clashes**<br>
    Different compilers follow different naming conventions. In
    particular Fortran compilers often (but not always) append or
    prepend one or two underscores to symbols in the object files. This
    is fine when resolving against other Fortran generated symbols, but
    requires coordination with the symbol names produced in C.

    A solution with an array of compile flags controlling the naming in
    the Fortran compiler, and underscores liberally scattered through C
    files is fragile and unreliable. It also makes it difficult to call
    library routines written in C.

-   **No checking of parameters**<br>
    The linker is extremely stupid - it only matches by parameter name.
    If this method is used, absolutely no checking is done on the
    parameters passed to the C routine from Fortran. This is a recipe
    for disaster, and will generate only runtime errors.

-   **Calling**<br>
    By default C passes arguments by value, whereas Fortran passes them
    by pointer. This requires writing wrappers for almost any
    non-trivial C library routine to access it from Fortran. Some
    constructs simply cannot be emulated.

-   **Variable types**<br>
    As an extension of the lack of checking of parameters, there is no
    checking of argument types across the Fortran/C interface. This
    relies on the Fortran and C code using the same types - in
    particular the same size of floating point and integer variables.
    This is extremely hard to guarantee, and these can fluctuate
    according to compiler flags. This can result in compiler and
    computer specific runtime errors that are extremely difficult to
    track.

All of these problems can be solved using structured interfacing, at the
cost of introducing a dependency on the Fortran 2003 standard.

All access to C routines in NECI must be through a declared interface.
This should be declared only once, in a module. An example is given
here:

```Fortran
interface
    ! Note that we can define the name used in fortran code, and the C
    ! symbol that is linked to independently.
    subroutine fortran_symbol(arg1, arg2) bind(c, name="c_symbol")
        use, intrinsic :: iso_c_binding, only: c_int, c_bool
        integer(c_int), intent(inout) :: arg1      ! Passed by pointer
        integer(c_bool), intent(in), value :: arg2 ! Passed by value
    end subroutine
end interface
```

A good summary of the rules and procedures for using this
interoperability are given in this Stack Overflow answer:
<http://stackoverflow.com/tags/fortran-iso-c-binding/info>

C++ routines can be made suitable for access from Fortran by prepending
symbol declarations in the C++ code with `extern "C"`.
