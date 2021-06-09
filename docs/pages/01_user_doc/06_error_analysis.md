---
title: Performing error analysis
---

## Performing error analysis

Data from an FCIQMC calculation is usually correlated. As a result,
standard error analysis for uncorrelated data cannot be used. Instead we
perform a so-called blocking analysis (JCP 91, 461). In this, data is
grouped into blocks of increasing size until the data in subsequent
blocks becomes uncorrelated, to a good approximation.

A blocking analysis can be performed in NECI in one of two ways.
Firstly, a rough blocking analysis is performed automatically after a
job is finished. The final result is output to standard output and
further information about the blocking analysis at various block sizes
is output to separate files, such as `Blocks_num` and `Blocks_denom`.
This should only be used as a rough and quick estimate as there are
issues with this approach. For example, the analysis starts as soon as
the shift is turned on. This is before the population has stabilised,
and so unusual results can occur in the analysis of the denominator and
numerator. Also, data is not taken from the optimal block size.

A better approach for a more careful analysis is to use the blocking
script in the utils directory, called blocking.py. The key command is

```bash
./blocking.py -f start_iter -d24 -d23 -o/ FCIMCStats
```

This will perform a blocking analysis starting from iteration
`start_iter`. The analysis should be started only once the energy
estimate, (column 11 in `FCIMCStats`) and the numerator and denominator
(columns 24 and 25) have stabilised and are fluctuating about some final
value. Just because the energy looks stable, it does not mean that the
populations is not still growing!

`-d24 -d23'` tells the script to perform the blocking on columns 25 and
24 of the `FCIMCStats` file, which correspond to the numerator and
denominator of the energy estimator, respectively. `-o/` tells the
script to also provide data for the results of dividing columns 25 and
24, which gives the energy estimate that we want.

Running this will produce a graph of the errors for both the numerator
and denominator as a function of the number of blocks (and therefore of
the block size). As the block size increases, the error estimates should
increase, tending towards the true values. Eventually the estimates will
plateau. This indicates that, at this block length, the data in the
blocks are uncorrelated to a good approximation, and the error estimate
calculated is accurate. The data from this block length should therefore
be used.

Each estimate of the error will also have an error on it. As the block
length increases this ‘error on the error’ will increase. One should
therefore use the *first* block length where the plateau is reached, so
as to minimise the error on the final error estimate.

If no plateau is seen in the plot then the simulation has not been run
for long enough, and needs to be continued by restarting from the
`POPSFILE`. It can take on the order of \(10^5-10^6\) iterations to
perform an accurate blocking analysis.

The `blocking.py` script will also output the final estimates on the
energy at the different block lengths. You should find the blocking
length where the errors plateau and read of the final estimates (the
rightmost columns) from here.

More information (including example plots, similar to those that
`blocking.py` produces) is available at JCP 91, 461.
