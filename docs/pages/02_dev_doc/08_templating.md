---
title: Code templating
---

## Code templating

NECI supports two ways of templating code, the python-based Fortran
preprocessor `fypp` (<https://github.com/aradi/fypp>) and the custom
script `tools/f90_template.py`. It is strongly suggested for new code to
make use of `fypp`, which allows for handling preprocessor flags using
python syntax, and provides an easy route to templating code with little
effort. All files named `*.fpp` have the `fypp` preprocessor applied. An
in-depth documentation can be found at
<https://fypp.readthedocs.io/en/stable/>.

As `fypp` support is a recent addition, most of NECI’s templates are
contained in `*.F90.template` files, which are automatically converted
by `tools/f90_template.py` into files with the corresponding names
`*.F90` prior to running the C preprocessor.

This mechanism exists to allow general code to be written, and reused
for different types, and combinations of types. While this is largely a
combinatorial pattern-matching and substitution problem, the templater
contains specific additional features to facilitate dealing with array
types in Fortran. There are also a number of specific considerations
that need to be made.

The templated code in NECI has been written largely by Simon Smart and
Alex Thom, who should be able to help with any particularly nasty issues
arising.

### How it works

Fortran permits multiple routines to be referenced by the same name
through the use of interface blocks such as

```Fortran
interface sub_name
    module procedure actual_name_1
    module procedure actual_name_2
end interface
```

which allows either of the routines `actual_name_1` or `actual_name_2`
to be called using the Fortran symbol `sub_name`. Note that these
procedure constructs can be used directly in the code for a hard-coded
set of routines which can be called from one interface name if desired
(see section
<a href="#sec:manualrenamingofroutines" data-reference-type="ref" data-reference="sec:manualrenamingofroutines">1.7.7</a>).

It is perfectly acceptable to have multiple interface blocks for a
specific routine name, so long as all of the referenced routines have
different calling signatures. That is, they must accept differently
typed arguments so that it is possible for the compiler to determine *at
compile time* which of the routines should be called. In principle the
actual routine names can always be used.

The Fortran templater creates one module per specified configuration,
each with a unique module name. It then performs substitutions on a
specified template model to create routines for all of the specified
combinations of input types. These routines have their names adjusted to
make them unique for each configuration, and an interface block is
created to make them accessible under their original name. Finally all
of these newly created modules are collected, with `use` statements,
into a macroscopic module which can be used from elsewhere.

### Overall structure

A sample templated module structure is given here for reference. The
different sections are explained below.
```Fortran

# This is the configuration block. Note that it has *.ini syntax,
# and that comments are preceeded by hashes.
[int]
type1=integer(int32)

[float]
type1=real(dp)

===================

#include "macros.h"

module module_name

    ! This is the module which is templated to generate the ensemble
    ! of routines with differing types
    use constants
    implicit none

contains

    elemental function test_fn(arg) result(ret)

        %(type1)s, intent(in) :: arg
        %(type2)s :: ret

    end function

end module


supermodule module_name
    !
    ! Here we include code that should be included in the module but
    ! does not need to be templated.
    !
end supermodule
```

### Configuration names and substitution

The top section of the file defines the templated configurations. It has
the structure of an INI file, and is processed by the standard python
ini file parser.

Configurations are defined by a name, contained in square brackets, and
then by a series of key-value pairs. All values are treated as strings
for the purpose of substitution in the main body of the routine.

Configurations are *inherited*. That is to say that all key-value pairs
(with the exception of `conditional_enable`) are carried forward to the
next configuration in the file unless they are overridden. This permits
quite sparse configuration files, at the expense of being a bit more
tricky to modify.

The length of configuration names must be considered. They will be
appended to module names and the subroutine and function names contained
therein. The templater does not have a magic means to circumvent the
Fortran 95 limit of 31 characters in any symbol. Therefore it makes
sense to use highly abbreviated configuration names.

The templater modifies the names of subroutines and functions as it
processes the module. As such, it is a little pick about syntax. Normal
routine decoration specifiers such as `pure` and `elemental` are
supported, but functions *must* be declared using the `result()`
specifier to define the return type.

If the special key `conditional_enable` is present, this is used to wrap
the generated module in `#if #endif` elements. See
`lib/quicksort.F90.template` for examples.

Values from the key-value pairs are directly substituted into the
templated module below, where they replace the element `%(key)s`. (This
specifier is the standard python named-string specifier). Keys named
beginning with `type` are treated specially, as described in the
following section.

The templater cannot circumvent Fortran line length limits. If necessary
a value to substitute can be extended over multiple lines by ensuring
the first character on a new line is a space, and then just continuing.
Make sure you remember the Fortran line continuation characters, as in
this example from the MPI wrapper code:

```Fortran
mpilen=((ubound(v,1)-lbound(v,1)+1)*(ubound(v,2)-lbound(v,2)+1)*&
         (ubound(v,3)-lbound(v,3)+1))
```

There is a special variable, which can be accessed using `%(name)s`.
This contains the name of the current configuration.

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

```Fortran
[int]                    =>   subroutine example_int(arr, elem)
type1 = integer(int32)   =>       integer(int32) :: arr(:)
                         =>       integer(int32) :: elem
                         =>       ...
                         =>   end subroutine
```

However, it may be necessary for the value which is being considered in
the array to itself be an array. An example of this would be the bit
representations used in NECI — a list of of these is a two dimensional
array, and any intermediate values would be arrays themselves.

In this case, an array type should be specified using the `dimension`
keyword, and the code will be automatically adjusted as follows:

```Fortran
[arr_int64]                            =>   subroutine example_int(arr, elem)
type1 = integer(int64), dimension(:)   =>       integer(int32) :: arr(:, :)
                                       =>       integer(int32) :: elem(:)
                                       =>       ...
                                       =>   end subroutine
```

Essentially the number of `:` delimiters appearing in the variable
definition is combined with the number of dimensions specified in the
type.

As a special case, temporary variables can be created of an appropriate
size which are either scalars, or have one dimension. For the definition

```Fortran
%(type1)s :: arr(:)
%(type1)s :: tmp(size(arr(1)))
```

then adjustment occurs as follows

```Fortran
[int]                    =>     integer(int32) :: arr(:)
type1 = integer(int32)   =>     integer(int32) :: elem
                         =>     integer(int32) :: tmp
```
And

```Fortran
[arr_int64]                           =>   integer(int32) :: arr(:, :)
type1 = integer(int64), dimension(:)  =>   integer(int32) :: elem(:)
                                      =>   integer(int32) :: tmp(size(arr(1)))
```


In a similar way, the references made to these variables within the
routines must be adjusted. This is to ensure that correct sized array
slices are used at all times. For the original templated code

```Fortran
    arr1(j) = arr2(i)
```

the following will result in the the templated output if the variable is
of an adjustable type and declared at the top of the function:

```Fortran
[int]                      =>    arr1(j) = arr2(i)
type1 = integer(int32)     =>
```
And

```Fortran
[arr_int64]                             =>     arr1(:, j) = arr2(:, i)
type1 = integer(int64), dimension(:)    =>
```

### The supermodule

In many modules, there are routines that do not need to be templated for
different variable types. As an example, within the MPI wrapper
routines, the code to initialise MPI is not variable type specific.

Code which is placed in the supermodule is not templated, and is
included directly in the final generated module.

### Optional parameters and lines of code

It is good practice to write templated routines as generally as
possible. This likely involves adding more functionality than is needed
in all cases, and switching this functionality on and off in some way.

For example, the sorting routine can sort multiple arrays in parallel,
according to the order in the first array (such as sorting a list of
determinants into energy order, where the energies are stored in a
separate array). It also needs to have comparison functions defined for
scalars as well as arrays.

The extent to which interesting features can be developed is limited
only by the developers imagination in using the template substition. But
two tricks are generally useful.

-   **Additional optional arguments**<br>
    Subroutines can easily be given flexible numbers of arguments. This
    is useful for adding additional functionality (and allows multiple
    templated routines to use the same `type` values). The templated
    subroutine definition
    ```Fortran
    subroutine example(arg%(extra_args)s)
    ```
    will generate the following code
    ```Fortran
    [simple]                    =>
    extra_args =                =>     subroutine example(arg)
                                =>
    [extended ]                 =>
    extra_args = , arg2, arg3   =>     subroutine example(arg, arg2, arg3)
    ```
    The next trick is useful for adding the type definitions of these
    additional arguments, and enabling the code which uses them.

-   **Switching off lines of code**<br>
    Lines of code in Fortran are trivially disabled when they are
    commented out. Prefixing lines with a switch-value allows it to be
    disabled. For example
    ```Fortran
    %(use_type2)%(type2) :: val()
    ```
    will allow an additional type to be used in a routine depending on
    the configuration:

    ```Fortran
    [unused]                        =>
    type2 =                         =>
    use_type2 =!                    =>        ! :: val ()
                                    =>
    [arr_real]                      =>
    type2 = real(dp), dimension(:)  =>        real(dp) :: val(:)
    use_type2 =                     =>
    ```

### Manual renaming of routines

The user can additionally manually create interface blocks for the
templated routines. This is useful where there is more than one possible
function to call for each of the variable types.

An example of this is given in the MPI wrapper functions, where there
are versions of routines that require manually specifying the lengths of
various parameters, and automatic versions which take the lengths from
the sizes of the arrays passed in. At the top of the templated module
definiton lie interfaces blocks such as

```Fortran
interface MPIReduce
    module procedure MPIReduce_len_%(name)s
    module procedure MPIReduce_auto_%(name)s
end interface
```

which makes use of the special `%(name)s` element to reference the
generated routines after templating.

In this case, the templated routines `MPIReduce_len` and
`MPIReduce_auto` will be available to the user as usual, but the
routines can both be called by the more generic name of `MPIReduce` with
the appropriate arguments supplied.

### Examples

All of the features of the templating code have been heavily used the
Shared Memory code, in `lib/allocate_shared.F90.template`, the sorting
code in `lib/quicksort.F90.template` and the MPI wrapper code in
`lib/Parallel.F90.template`. Other less aggressively used case can be
found elsewhere.
