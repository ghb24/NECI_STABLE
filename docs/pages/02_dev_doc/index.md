---
title: Developer's guide to NECI
---

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
