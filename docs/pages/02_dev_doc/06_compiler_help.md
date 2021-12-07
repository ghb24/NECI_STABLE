---
title: Let the compiler help you
---

## Let the compiler help you

If real NECI is compiled in debug mode, warnings will be treated as
errors. This means that code with real numbers may not produce any
warnings to be merged into main development branches. The complex
version of the code however, does not treat warnings as errors.

Unfortunately the warnings for unused variables had to be deactivated,
because there are too many incidents.

@warning
There are too many conversion warnings in the complex NECI code that
could not be cleaned up yet and probably lead to serious bugs. It is
necessary to clean them first before complex NECI can be used reliably.
@endwarning

Sometimes a warning is a false-positive. To work around such problems
there is a `WARNING_WORKAROUND_` compile flag that gets activated, if
warnings are activated.

A common false-positive warning is an unused variable that has to stay
in the code, because it is e.g. in the interface of a function that is
the target of a function pointer. For this case the `unused` macro
exists. It is not necessary to put this macro behind the
`WARNING_WORKAROUND_` compile flag. It is recommended to mark unused
variables directly after declaration to make it explicit to the human
reader.

```Fortran
#include "macros.h"
integer, intent(in) :: arr1(:), n
real(dp), inten(in) :: arr2(:, :)
integer :: ierr

unused(arr2); unused(n)
```
