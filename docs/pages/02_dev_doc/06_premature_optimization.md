---
title: Don’t optimise prematurely
---

## Don’t optimise prematurely

Obviously, good algorithm design is important. If an algorithm scales
badly, then no implementation will be able to salvage it.

However, there are many tricks and optimisation that can be made to eek
out small and large performance gains in the implementation of a
particular algorithm. It is important not to optimise too early for a
number of reasons.

-   Good optimisation is extremely time intensive. On the whole time is
    better spent getting the code to work, and making the algorithm
    efficient. Once an implementation works, then code can be profiled
    and performance improved.

-   Optimisations are often highly non-obvious, involving storing
    information in unexpected ways and places, leading to code that is
    hard to write, harder to read and extremely bug prone.

-   The compiler is very good. The obvious ‘tricks’ that you see will be
    done by the compiler anyway.

-   One place where performance gains can legitimately be made is in
    avoiding conditional switching. In many cases this would involve
    duplicating code paths, and horrifically breaking the DRY principle
    above. There are occasions that this is worthwhile, but this should
    be actively justified by profiling data rather than just a hunch.
