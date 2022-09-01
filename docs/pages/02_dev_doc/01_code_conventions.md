---
title: Code conventions
---

## Code conventions

The code in NECI has been developed over a number of years by many
different developers, and has very little standardisation of approach or
code appearance. This is not an example to copy!

We are trying to (gradually) normalise sections of code, and isolate
those sections which are old and generally unredemable from the rest of
the code base. As such there are a number of restrictions we place on
code in NECI, and a range of other guidelines.

**Fortran standard**<br>
Due to the use of procedure pointers, a reasonably up to date compiler
supporting (at least some of) the Fortran 2003 standard is *required* to
compile NECI. The C-interoperability and procedure pointer features of
Fortran 2003 should be used. Other features of this standard should be
used sparingly, as Fortran 2003 support in compilers is patchy at best.

Otherwise, code should be written to the Fortran 90/95 standard. In
particular, several features of FORTRAN 77 should be avoided at all
costs:

`DO` statements using `REAL` type loop variables.

Assigned `GOTO` statements

Cray pointers (declared with the format `pointer (ptr, pointee)`)

Implicit variables. `implicit none` MUST appear in every module, or
interface statement.

Implicitly typed routines and subroutines. See section on modules and
interfaces.

`COMMON` blocks for sharing data between files.

All of these features do, or have appeared, in NECI at some point. It is
also highly advised that `intent` arguments should be used for all
argument declarations. This both improves performance, and the ability
of others to quickly identify what the routine is attempting to do.

Regarding code layout, the PEP8 guidelines of the python language lead to
well readable code and are mostly applicable to Fortran as well.
(<https://www.python.org/dev/peps/pep-0008/>)

**CAPITAL letters**<br>
Fortran is a case insensitive programming language.

For historical reasons a large proportion of FORTRAN 77 code was written
entirely in CAPITAL LETTERS (with the exception of displayable strings).
This is extremely bad practice.

Humans generally read by recognising word shape. This is obliterated in
fully capitalised text, making code much harder to read, and typos
especially difficult to identify.

**Indentation**<br>
Indentation of sections of code should use spaces (and not tabs). The
Fortran 95 standard explicitly rejects the use of tabs, and tabs in
source code will elicit warnings from the compiler.

All indentations should be multiples of 4 characters.

Source code in .F files (old-style FORTRAN 77) has specific layout
restrictions. In particular an initial indent of 7 spaces. This style
should not be mimicked elsewhere.

**Code line length**<br>
The Fortran 90/95 standard restricts line lengths to a (hard) maximum of
132 characters. Code with lines longer than this *may* work on *some*
compilers, but this limit should be avoided.

This limit applies after preprocessing has been applied. A number of our
macros in `macros.h` can create lines of considerably longer length if
not used carefully. These may require using temporary variables with
shorter names to control the line length.

The fixed-format FORTRAN 77 code is restricted to 72 characters per
line.

On a 19" monitor at standard resolution, two columns of code vertically
split and side by side use 79 characters each. This is a convenient
soft-limit to use - although it is not trivially achievable in all code,
and overall readability should be prioritised.

**Variable name conventions**<br>
There are a number of competing conventions for variable and function
names within NECI. That said, there are a number of existing conventions
that it is *useful* to be aware of and which new code should keep in
mind.

`CamelCase` or `snake_case` should be used to provide descriptive
variable names. Prefer `snake_case` when reasonable. The wider the scope
of a variable, the longer and more descriptive the name should be.
Trivial local variables (loop indices, etc.) can and should be trivially
named. Variables should not, ever, be used entirely in capital letters.

`tVariableName` is a logical control variable. Normally globally
declared in a module for switching on (or off) an overall feature, or
signalling overall calculation state.

`TypeName_t` is a user-defined type.

`nVariableName` is an integer containing a count of a number of a given
entity.

`Variable`, `AllVariable` are paired sets of variables tracking an
extensive property of a simulation. That is properties which can be
accumulated on an individual node, but that the system-relevant property
needs to be collected from all nodes and amalgamated. This is done once
per iteration or once per update cycle as appropriate.

`CamelCase_t` CamelCase and a trailing t denote a derived type.

Fortran 95 restricts variable names to 31 characters. Although Fortran
2003 extends this to 63, making use of this extension can cause problems
with some compilers, and this should be avoided.

**Subroutine decoration (especially intent statements)**<br>
Subroutine and function declarations should be decorated to the greatest
extent feasible. This should restrict the variables to only their
expected role in a function.

In particular, all function arguments should be decorated with either
`intent(in)`, `intent(out)`, or `intent(inout)` as appropriate. When
absolutely necessary, use `value` to pass arguments by-value. The
additional decorations `optional` and `target` can be used with care.

The supplied arguments should be as restrictive as possible, to maximise
the likelihood of the compiler catching programming errors.

The single exception to this rule is for routines that are to be stored
in procedure pointers. These routines must exactly match the definition
of the relevant `abstract interface`, which may be more general than is
required for the specific case.

The arguments should be sorted by non optional `in`, `inout`, `out` and
then optional arguments. The declaration of dummy arguments should
appear in the same order as the argument list. There should be an empty
line between dummy argument declarations and local variable
declarations.

If a procedure is often, but not always, called with the same argument
think about making it `optional` using the `def_default` macro and
introduce a placeholder local variable with appended underscore.

If possible add the `pure` attribute. This shows the human that there
are no side effects and makes parallelization and encapsulation easier.

If a function is `pure` and operates on scalars, add the `elemental`
attribute to automatically map it elementwise onto arrays.

The following toy function is pure, can be applied onto arrays and
scalars alike. If the exponent is ommited, it defaults to squaring.

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

**Data types**<br>
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

**Array declarations**<br>
Fortran arrays can be declared in multiple ways. In particular, the
dimensionality of an array can be declared on the variable itself, or as
part of the type declaration;

```Fortran
integer, dimension(10, 20) :: arr1
integer :: arr2(10, 20)
```

In general the latter declaration is preferred for two reasons:

1.  It is clear that the array property is attached to the variable, and
    not to the type. When scanning data declarations it is not possible
    to mistake a scalar for an array.

2.  Multiple different arrays can be declared in the same data
    declaration with different bounds.

The only exception to this is in templated code, where the bounds of
arrays need to be varied.

Where arrays are passed as arguments to a routine, they can be passed in
three ways

```Fortran
integer, intent(inout) :: arr(*)
integer, intent(inout) :: arr(10)
integer, intent(inout) :: arr(:)
```

The first form should be avoided wherever practical, as it prevents any
knowledge of the array dimensions being carried into the code. This
means that whole-array manipulations will no longer work.

The second two types can be used in different circumstances. The first
essentially overrides the array dimensions passed in. This can be useful
for re-indexing arrays (e.g. treating a zero-based array as one-based).

The last approach allows a receiving routine to inspect the array bounds
as passed in by the calling routine. This maximises the extent to which
the compiler and debugging tools can assist in finding errors in the
code, and should be used wherever possible.

**use statements**<br>
Globally declared symbols can be shared between modules using `use`
statements. Generally, specific symbols should be included rather than
all symbols in a module using the notation
`use module_name, only: symbol, ...`.

Modules containing only data that have been separated for the purposes
of dependency resolution can be fully included into their related
modules (e.g. `CalcData` and `Calc`).

Use statements should (where possible) be located in the module header,
and not in individual subroutines - this avoids some serious issues
associated with compilers resolving conflicting dependencies. If the
same thing is included in multiple places in a file, the compilers
dependency resolution tree can become very large, and use a lot of time
and memory to resolve unambiguously.

**ASSERT statements**<br>
`ASSERT` is a macro, defined in `macros.h`. In an optimised build these
statements are entirely removed, and in a debug build they will cause
execution to be aborted with an error message if the condition specified
is not met.

In NECI, the error message contains the current file and line number. It
also includes the current function, which must be manually supplied in a
constant named `this_routine`.

An example assert statement, in a function that takes an array with the
same number of elements as there are basis functions, would be:

```Fortran
subroutine foo(arr)
    integer, intent(inout) :: arr(:)
    character(*), parameter :: this_routine = 'foo'
    ASSERT(size(arr) == nBasis)
    ...
end subroutine
```

@warning
Be careful not to use tests with side effects in `ASSERT` statements. In
the optimised build the tests will not be called, and this can introduce
bugs that appear in only one of the optimised or debug builds.
@endwarning

**Floating point comparison and integer division**<br>

It is usually a bad idea to test floating point numbers for
(in-)equality using `==` or `/=`. Equality should rather be tested with
an expression like \(|a - b| < \epsilon\). In the `util_mod` module
there are `near_zero` and `operator(.isclose.)` which should be used for
this purpose.

If one divides two integers `5 / 3 == 1` the result gets truncated to
the nearest integer. Sometimes this is not wanted, and the compiler
warns about it. For this reason one should use `5 .div. 3` to make it
explicit that integer division is indeed wanted.

**Tools for adhering to the style guide**<br>

@warning
It is better to write nice code from the beginning on instead of relying
on automated tools. This section is meant for existing code, that has no
consistent indentation and other flaws.
@endwarning

One recommended program to prettify Fortran free-format code is
`fprettify`. It can be installed with `pip3 install fprettify` and is
automatically installed on the Alavi workstations. The `--user` option
might be required for installation if you do not have `sudo`-rights.

The NECI codebase already contains the correct configuration files, so
it is sufficient to just call `fprettify` on a file in `src/` or
`src/lib`.

**Operator Layout**<br>
The following guidelines are recommendations for formatting of code and
represent the configuration of the `fprettify` tool explained in the
previous paragraph. These binary operators should be surrounded with a
single space on either side: assignment (`=`), comparisons
(`==, <, >, /=, <=, >=`), Booleans (`.and., .or., .not.`).

If operators with different priorities are used, consider adding
whitespace around the operators with the lowest priority(ies) to make
the expression readable easily. Use your own judgment; however, never
use more than one space, and always have the same amount of whitespace
on both sides of a binary operator.

```Fortran
! Recommended
    i = i + 1
    x = x*2 - 1
    hypot2 = x*x + y*y
    c = (a+b) * (a-b)

! Also possible
    x = x * 2 - 1
    hypot2 = x * x + y * y
    c = (a + b) * (a - b)

! Not recommended
    i=i+1
```

Line breaks should happen before binary operators for easy association
of operator and operand

```Fortran

! Not recommended: operators sit far away from their operands
    income = gross_wages + &
             taxable_interest + &
             (dividends - qualified_dividends) - &
             ira_deduction - &
             student_loan_interest

! Recommended: easy to match operators with operands
    income = gross_wages &
             + taxable_interest &
             + (dividends - qualified_dividends) &
             - ira_deduction &
             - student_loan_interest
```

Please use the new C-style relational operators.
```Fortran
! Recommended
      ==   /=   <    <=    >   >=
! Not recommended
      .EQ. .NE. .LT. .LE. .GT. .GE.
```

**Whitespace in Expressions**<br>
Avoid extraneous whitespace in the following situations.
```Fortran

! No whitespace immediately inside parentheses:
    Yes: spam(ham(1), f(eggs, 2))
    No:  spam( ham( 1 ), f( eggs, 2 ) )
! No whitespace immediately before the open parenthesis that starts
! the argument list of a function call or array indexing:
    Yes: spam(1)
    No:  spam (1)

! No whitespace immediately before a comma, semicolon, or colon:
    ! Yes:
        use module, only: cool_function
        integer, allocatable :: A(:, :)
        integer, allocatable :: A(:,:)
    ! No:
        use module , only : cool_function
        integer , allocatable :: A(: , :)
! No whitespace around the = sign when used to
! call a function with a keyword argument.
    ! Yes:
        pow(2, n=3)
    ! No:
        pow(2, n = 3)
```

In a slice the colon acts like a binary operator, and should have equal
amounts on either side (treating it as the operator with the lowest
priority). In an extended slice, both colons must have the same amount
of spacing applied. Exception: when a slice parameter is omitted, the
space is omitted.

```Fortran
! Recommended
    ham(1:9), ham(1:9:3), ham(:9:3), ham(1::3), ham(1:9:)
    ham(lower:upper), ham(lower:upper:), ham(lower::step)
    ham(lower+offset : upper+offset)
    ham(: upper_fn(x) : step_fn(x)), ham(:: step_fn(x))
    ham(lower + offset : upper + offset)

! Not recommended
    ham(lower + offset:upper + offset)
    ham(1: 9), ham(1 :9), ham(1:9 :3)
    ham(lower : : upper)
    ham( : upper)
```




**Contained procedures**<br>
From Fortran2003 onwards it is possible to define procedures inside
procedures. The inner procedure has access to the local scope of the
outher procedure. (Similar to closures in other languages.) These
contained procedures allow to cleanly eliminate some reasons, why one
would like to use a global variable. Besides it also leads to less
wrapper functions that are similar, but not the same. This can be used
both to avoid code duplication and to avoid passing unnecessary
arguments.

Let’s assume there is a `fancy_function` with ten arguments. One of them
is a `real, intent(in) :: x`. If `fancy_function` is called several
times with only a varying `x` there are many lines of code doing more or
less the same. In former versions of Fortran this usually lead to the
introduction of a wrapper function `my_fancy_function` that depends
explicitly only on `x` and gets the nine other arguments using global
variables that have to be defined before calling `my_fancy_function`.

Another code that uses `fancy_function` keeps `x` constant, but varies
another argument `y`. This leads to a second wrapper
`my_fancy_function_2` and more global variables. In addition both
wrapper function might be used only at one place.

If one instead defines the wrapper function as internal procedure only
where it is used there is no need for a wrapper function. This resembles
the process of *currying* in functional languages.

```Fortran
function calculate_something(a) result(res)
    real, intent(in) :: a
    real :: res

    ! arg2 to arg10 are visible identifiers here

    res = exp(f(a)) + 3

contains

    function f(x) result(res)
        real, intent(in) :: x
        real :: res

        res = fancy_function(x, arg2, arg3, ..., arg10)
    end function

end function
```

Another use case for contained procedures is to extract recurrent parts
of a subroutine/function without exposing them to module scope.

A contained procedure cannot make use of the `contains` statement.

**Modules and interfaces**<br>
It is an aim to make the dependency between different code parts as
small, as unidirectional and as explicit as possible. Module are a great
tool to achieve that goal.

A good example are the different excitation generators of NECI. It is
possible to seamlessly exchange different excitation generators because
they all have the same public interface. On the other hand they greatly
differ in their implementation details and each excitation generator
uses different helper functions. It is possible to use an excitation
generator without knowing about the implementation details and helper
functions. It is even advised to not rely on or assume any
implementation detail for a specific excitation generator.

There are some rule of thumbs to achieve similar results in the
architecture of other code. If one is implementing a new module it is
good to start with the `private` keyword to make all identifiers private
to that module and explicitly thinking about which identifiers should be
accessable from outside and declare them `public`. These `public`
identifiers should only change for a good reason afterwards and should
be well documented. Variables that should have read-only access from the
outside (a computed energy for example) can be declared
`public, protected`.

If other modules are imported with `use, only:` then it is easy to see
on which code a module relies.

If there are generic functions it is possible to declare only the name
of the generic interface as `public` and keep the concrete
implementations `private`.

The use of `private` helper functions has the same benefit as contained
procedures. It allows to write ad-hoc wrappers that have access to the
module scope, but do not require `public` global variables.

If possible functions should be declared `pure` or `elemental`.

**Example module layout**<br>
A sample module layout is given below:

```Fortran
#include "macros.h" ! This enables use of our precompiler macros.
module module_name

    ! To the extent possible, include statements should be at the
    ! beginning of a module, and not elsewhere.
    ! If possible they should import only actually used identifiers.
    use SystemData, only: nel, tHPHF
    use module_data
    use constants
    implicit none

    ! To the extent possible, declare all identifiers of
    !   a module as private by default
    !   and export explicitly with the public keyword.
    private
    public :: sub_name, fn_name, calculated_energy
    ! Cannot be changed from the outside.
    protected :: calculated_energy

    real(dp) :: calculated_energy


    ! Add an interface to an external (non-modularised) function
    interface external_fn
        function splat_it(in_val) result(ret_val) &
                                  bind(c, name='symbol_name')
            ! n.b. interface statements shield from modular includes
            import :: dp
            implicit none
            integer, intent(in) :: in_val
            real(dp) :: ret_val
        end function
    end interface

contains

    [pure|elemental] subroutine sub_name(in_val, out_val)

        ! This is a description of what the subroutine does

        integer, intent(in) :: in_val
        real(dp), intent(out) :: out_val

    end subroutine [sub_name]


    [pure|elemental] function fn_name(in_val) result(ret_val)

        ! This is a description of what the function does

        integer, intent(in) :: in_val
        real(dp) :: ret_val

    end function [fn_name]

end module
```

**Error handling**<br>

It is very important to always be in a well defined state and to be
deterministic up to the stochastic noise of the Monte-Carlo simulation.
For this reason, errors that cannot be handled inside a procedure and
are not handled by calling code should crash the program. This
philosophy follows the convention of the fortran intrinsic `allocate` or
the `mpi_f08` subroutines. The error code should be optional and calling
code should only pass it, if they handle all possible error cases. The
called function should not abort the calculation, if an error code is
present and should abort the calculation if this is not the case.

Calling code should only ask for the error code if the return value is
handled as early as possible. Do not write

```Fortran
allocate(A, stat=ierr)
allocate(B, stat=ierr)
allocate(C, stat=ierr)
if (ierr /= 0) call stop_all(...)
```

but either check after each allocation (to know which array failed) or
just ommit `stat=ierr`. The intrinsic allocations stops the program if
it cannot allocate.


**One line if statements**<br>

One-line if statements should fit in one line.
If they have to be continued via `&`, please use a proper `then` ... `end if`
pairing.

```Fortran
! This is allowed
if (cond) statement


! This is forbidden
if (cond) &
    statement
```

The underlying reason is that the danger is too high, to go from
```Fortran
if (cond) &
    statement_1
```
to
```Fortran
if (cond) &
    statement_1
    statement_2
```
where people might overlook that `statement_2` is always executed.
