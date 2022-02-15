---
title: Debugging tips
---

## Debugging tips

Code breaks. Sometimes it appears to rot with time. Finding bugs takes
the majority of most programmers time, and practices which make this
quicker are invaluable!

Obviously, the techniques you will use will depend on the nature of the
problem (tracking down small numerical changes in output is generally
much harder than finding what causes a segfault), but there are number
of tricks which can help!

### Build configurations

As soon as you have a problem, build a debug rather than an optimised
version of the code. This cause a large array of changes to the compiled
code:

-   **Array bounds checking**<br>
    All Fortran arrays have well defined bounds on all of their
    dimensions. In debug mode the compiler will insert code to check
    that all memory accesses are within these bounds. If not, execution
    will be terminated with a message indicating at what line of what
    file the error occurred. If running in a debugger (see later)
    execution will be interrupted at this point.

-   **Disable optimisations**<br>
    Hopefully this will not make the bug go away! If it does, you are
    almost certainly looking at either an uninitialised variable, or
    access beyond the end of an array.

    The primary purpose of disabling optimisations is to make the
    mapping between the source code and the executable more linear. This
    results in any error messages, and the output of any tools, being
    easier to interpret.

-   **Adds debugging symbols**<br>
    When your code crashes it is really useful to know what routine was
    running, and what the stack trace (list of routines that have been
    called to get to this point in the code) is. Adding debugging
    symbols provides the information to convert the memory addresses
    into files and lines of source code. This makes error messages
    useful.

-   **Enables the ASSERT macro**<br>
    There are many consistency checks internally in NECI that can be
    turned on in debug mode. Particularly for difficult-to-find bugs,
    these are likely to fail substantially earlier in a run than it is
    possible to view the problems in the normal output.

-   **Defines \_\_DEBUG**<br>
    Any blocks contained inside `#ifdef DEBUG>` sections are only
    enabled in debug mode. Most of these contain either additional
    output specifically targetted to make debugging easier, or
    additional consistency checks.

### ASSERT statements are your friend

Generally the time to add `ASSERT` statements to your code is when you
are writing it. Debugging will make you acutely aware of their benefit.
If you have a suspicion in which bits of code a bug might lie, liberally
sprinkling it with `ASSERT` statements can help to find bugs.

More usefully, when you find the bug, add `ASSERT` statements to the
code to catch similar errors in the future.

### Learn to use your tools

Most of the programming tools in existence are for the purposes of
debugging. Learn to use them! Practice using them. The really work.

-   **Debuggers**<br>
    The debugger is the most powerful tool you have. Essentially you run
    your code in a harness, with the debugger hooked into everything
    important. On most Linux systems, the most readily available
    debugger is `gdb`.

    The debugger can trap any execution errors, and will interrupt
    (break) execution of your program at this point — while preserving
    all of the memory and execution state. You can examine the call
    stack (the nested list of functions that have been called), and the
    status of all of the memory. You can also change the contents of any
    relevant memory and continue execution to see the effects.

    Break points can be added to the code at any line of any function,
    and execution interrupted at that point. If your bug is only showing
    late in execution, you can interrupt on the \(n\)-th time a function
    is called. You can step through the execution line by line in the
    source code, examining all variables, and watch precisely what goes
    wrong.

    The debugger is the swiss-army sledgehammer of tools.

-   **Valgrind**<br>
    Valgrind is a tool for memory debugging. In essence it replaces most
    of the memory manipulation primitives provided by the operating
    system with instrumented versions. It will track what happens to
    memory, where it is created, where it is destroyed, and what code
    (in)correctly accesses it.

-   **Intel Inspector XE**<br>
    If you have access to Intel tools, the Intel Inspector is an
    extremely powerful debugger, memory analysis tool, performance
    enhancement and problem tracking tool. An top of a good interface to
    normal debugging tools, it can apply a lot of analysis to your
    code’s execution, and then filter through it to find where things go
    wrong.

### Machete Debugging

When bugs are particularly non-obvious, a good first step is to reduce
the complexity of the problem. Try and remove as much code as possible.
A good rule of thumb is to remove half of the code, and retest. If the
bug has gone away, then it required interacting with the other half of
the code.

This technique can very quickly isolate a test case (that still fails)
into something of manageable size. In doing so, the bug normally becomes
obvious. You can then revert to the original code, and fix the problem
there.

In a complex code, this process can be tricky. Lots of the sections of
code are coupled to each other. This can require some creative thinking.
Choose the domain of your problem carefully (it may only be necessary to
apply the machete to one file), and aggressively decouple sections by
commenting out function calls.

Remember — it doesn’t matter if you break functionality doing this (you
obviously will) so long as you retain the buggy behaviour in the code
that is left.

This is really good for finding strange memory interactions — left with
the two routines that are trampling on each others toes.

### fcimcdebug 5

For bugs that change the trajectory of simulations, turning on the
maximum debug output level in NECI (by adding `fcimcdebug 5` to the
LOGGING block of the input file in a debug build) can often give a quick
insight into where the problem is located.

Many problems will exhibit with obviously pathological behaviour. For
example trying to spawn `NaN` particles, or an extremely large number,
generally indicates a problem in either Hamiltonian matrix element
generation or calculating the generation probabilities. Similarly
walkers being spawned with repeated or zeroed orbitals are a give away.

Similarly, if the bug is recently introduced and a prior version
functioned correctly (such as a test code failure)

### Look at the git logs

If you are chasing down a regression (such as a testcode failure), then
there is certainly a version that used to work, and a version which now
does not. Frequently there are not many commits between these two
versions — have a look at what they are. Often the failure is really
obvious!

If there are a large number of commits between the last known good
commit and the first known bad commit, then using `git bisect` is a very
efficient way to locate the first bad commit and the last good commit.
