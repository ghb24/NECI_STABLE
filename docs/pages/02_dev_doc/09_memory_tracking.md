---
title: Tracking memory usage
---

## Tracking memory usage

NECI contains automated tracking of memory usage. This enables output
statistics to indicate which memory uses are dominating during a
calculation.

The `MemoryManager` module keeps track of all units of memory, and
assigns a tag value to each of them. It is the responsibility of the
developer to store this tag, and pass it when the memory is deallocated.
This tag is an `integer`.

Memory should always be allocated using error checking. That is, an
allocate statement should always be passed an error value as follows

```Fortran
integer, allocatable :: arr1(:)
real(dp), allocatable :: arr2(:,:)
integer :: ierr
allocate(arr1(10), arr2(20, 30), stat=ierr)
```

This value will be zero if the allocation was successful, and non-zero
otherwise. The memory logging routines check this value, and report an
error if the memory allocation failed.

Memory is logged using the functions `LogMemAlloc` and `LogMemDealloc`.
