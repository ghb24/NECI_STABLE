---
title: Quickstart Tutorial
---

## Basic NECI Tutorial

FCIQMC is not a blackbox method and as such may be daunting to first approach.
In addition, there are many variants which excel in different kinds of problems.
NECI has tools to make working with FCIQMC an easier task.
While NECI is written predominantly in Fortran, in order to use it you need not know any programming language.
In this tutorial, we assume an at least elementary understanding of the FCIQMC algorithm.[@Booth2009]

The goal of this tutorial is to provide a practical, brief supplement to the full [NECI user's guide](index.html), which serves us a much more detailed reference.
We will use NECI to calculate the ground-state energy of the Nitrogen dimer in a (6,6) active space.
We use the equilibrium geometry, 2.074 bohr radii.

@note
For this problem, FCIQMC will actually underperform compared to conventional Davidson CI.
In addition, Davidson CI has no stochastic error (unlike FCIQMC), so it is strongly preferred for this problem.
However, this is a relatively cheap example with which to get comfortable with NECI and FCIQMC.
FCIQMC has better scaling behaviour, can be better parallelised, and benefits from sparsity.
If you go to (e.g.) a (20, 20) CAS then FCIQMC outshines Davidson, which is not a feasible method for problems of this size.
@endnote

First, we must generate the FCIDUMP file which contains the information about 1- and 2-electron integrals.
NECI is a solver for the CI-problem and *not* a standalone quantum chemistry suite, and cannot do this itself.
For this, choose any program that can generate these (e.g. PySCF, Molpro, Molcas).
This has been done for you, and you may download the file [here](|media|/N2_neci_files/FCIDUMP).

@todo
Perhaps also include the Molpro code to generate this FCIDUMP file?
@endtodo

### Anatomy of an Input File

In order to run a NECI calculation, we must create an input file. Here is an example for the FCIDUMP file provided, called `n2_neci.inp`.
```
# comments are given like this
( or like this (which is deprecated)

# simple N2 FCIQMC calculation
# for more complex FCIQMC variations, see the keywords for additional options
# such as the Hubbard model, transcorrelated options, or GAS-CI

title

  # read integrals from FCIDUMP
  system read
    electrons 6
    nonuniformrandexcits pchb
  endsys

  calc
    nmcyc 10000
    # for reproducibility
    seed 8

    totalWalkers 50000
    tau 0.01 search

    # use the initiator method
    truncinitiator
    addtoinitiator 3

    methods
      method vertex fcimc
    endmethods

  endcalc

  logging
    hdf5-pops
  endlog

end
```

@note
these keywords are case insensitive.
@endnote

All these keywords (and plenty more) are explained in the [next section](03_calculation_inputs.html). Let's break down the structure of this input file.

- Comments can be added in the code with `#`. (A deprecated comment symbol found in legacy inputs is `(`)
- First, the actual input starts with `title`, which is mandatory, and must also end with `end` (i.e. wrap the program in this block).
- Next, we have the `system` block, which is also mandatory.
    - The `system` keyword has a mandatory argument which comes directly after it on the same line. Here, we use `read` (as in `system read`), as we are doing the FCIQMC calculation from an FCIDUMP file.
    - We must also specify the number of electrons in the system, with the `electrons` keyword, followed by the number of electrons. Since we are doing a (6,6) calculation, we have `electrons 6`.
    - We must terminate the `system` block with the `endsys` keyword so that NECI knows to stop looking for `system` keywords.
- Then, we have the mandatory `calc` block, which is necessarily terminated with `endcalc`. This block in particular has many options and potential keywords; here we use only a small subset.
    - We specify the number of iterations the FCIQMC calculation will do, with `nmcyc`. Specifically, we specify 10000 iterations.
    - We also specify a seed with the `seed` keyword. FCIQMC is a randomised algorithm, and setting the seed simply ensures we get the same result every time (useful for testing or checking stochastic effects, for example). This keyword is optional.
    - We also must include the `totalWalkers` keyword, followed by the target number of walkers. Once this number is reached, NECI will enter variable-shift mode; that is, the shift will vary so as to keep the number of walkers constant. We wish to have a statistically significant number of iterations where the number of walkers is roughly constant, as we will see. In this case, we have `totalWalkers 50000`.
    - We must also include `tau` which is the size of the time step per iteration. The additional keyword `search` is optional but useful for stability.
    - Another form of FCIQMC is i-FCIQMC, which uses "initiator" walkers to speed up the calculation.[@Cleland2010] This is optional, but generally recommended. Presence of the `truncinitiator` keyword in the `calc` block indicates i-FCIQMC, and `addToInitiator 3` means that any determinant with a population >= 3 will become an initiator.
    - Finally, as a mandatory subblock *inside* the calc block, we have the `methods` subblock, which necessarily ends with `endmethods`. Here, we simply specify what kind of calculation to do. We choose `method vertex fcimc`, which simply means to run an FCIQMC calculation.
- Finally, we have an *optional* `logging` block which ends with `endlog`. By default, NECI will keep track of the population of walkers in a "POPSFILE", which is by default ASCII. Here, we wish to have an HDF5 POPSFILE which is generally better-performing. To do this, inside the `logging` block we have the `hdf5-pops` keyword.

@note
You cannot use the `hdf5-pops` keyword if you did not build NECI with HDF5.
@endnote

### Running NECI

After building, the NECI executable will be in `path/to/neci/build/bin/neci` (e.g. if I installed NECI in my home directory, it would be `~/neci/build/bin/neci`).

NECI can be run directly as any executable:
```bash
path/to/neci/build/bin/neci n2_neci.inp
```
but parallel execution is usually desired.
To run the above input file in parallel, we must use the respective MPI commands (`mpirun`, `mpiexec`, etc.)
```bash
mpirun -np 4 path/to/neci/build/bin/neci n2_neci.inp
```
where you can replace the `4` with however many processors with which you want to run (4 being a *very* modest number).
This will print a lot to standard output, which you may wish to capture, e.g.
```bash
mpirun -np 4 ~/neci/build/bin/neci n2_neci.inp > n2_neci.out
```

@note
If you made a mistake somewhere in your file (for example, a typo), NECI tries to give useful error messages.
However, MPI can sometimes interferes with this. Before doing a full calculation, you may want to run simply without `mpirun`, e.g. simply
```bash
~/neci/build/bin/neci n2_neci.inp > n2_neci.out
```
If this starts to run without an error, then you may stop it and run with `mpirun`.
@endnote

@note
NECI has a useful utility to dynamically change variables *while it is running*.
To do this, create a file called `CHANGEVARS`, and input whatever keyword you wish to change,
  e.g. if you think `nmcyc 10000` is too small, enter into `CHANGEVARS` the text `nmcyc 20000` (or whatever else).
  Then, the value will be updated in the next iteration of NECI. This is helpful for interacting with a running simulation.
@endnote

@warning
It is generally not advisable to terminate a program with CTRL+C. Instead, use the above `CHANGEVARS` trick, to create a soft exit. Simply open/create `CHANGEVARS`, type in `SOFTEXIT` and save. Alternatively, you can do this in one line on the terminal with
```bash
echo SOFTEXIT > CHANGEVARS
```
@endwarning

This calculation should produce a few other files in the directory, namely:

- `popsfile.h5` (or `POPSFILE` if you did not include `hdf5-pops`)
- `INITATORStats`
- `FCIMCStats`

The popsfile will be very useful in case we wish to continue running the simulation. `FCIMCStats` has columns of useful data, which we will explore now.

### Checking Convergence and Analysing Results

The file `FCIMCStats` has several useful columns which you will want to plot to ensure convergence.
To do this in one line, [there is a convenience script here](|media|/N2_neci_files/plot_fcimcstats.plt),
which is run with gnuplot via `gnuplot plot_fcimcstats.plt` and will output plots of the most useful columns to a new `plots/` directory.
Your results should look something like this:

1. ![](|media|/plots1/check_totE.png)
2. ![](|media|/plots1/check_totWalkers.png)
3. ![](|media|/plots1/check_refPop.png)
4. ![](|media|/plots1/check_shift_energy.png)
5. ![](|media|/plots1/check_numerator.png)
6. ![](|media|/plots1/check_denominator.png)

Plot (1) is the most immediately useful plot, as it gives you a quick estimate of the total energy from the calculation (namely, you can see we have around 108.98 Hartree). For all of these, we expect variable behaviour until the total walkers (plot (2)) reaches the target total walkers (as per our input file, 50000 in this case). Then, we want all six of these plots to roughly plateau, as they all do above. Furthermore, FCIQMC has a built-in consistency check whereby the energy can be calculated in two independent ways: namely, by the shift (which is only updated once the target number of walkers is reached) and the projected energy. These are plotted on top of each other in plot (4). As we see, they agree once we have run for a long enough time.

Once we are confident that all these plots exhibit plateaus for sufficiently large step numbers, we proceed with an error analysis. However, since FCIQMC calculations generally have correlated data, we cannot use standard error analysis, and here we use blocking analysis.[@Flyvbjerg1989] A script to do blocking analysis is included in the NECI repository: `path/to/neci/utils/blocking.py`.

@todo
The blocking script is still written in Python2, and IMO should be updated to Python3 (which should be easy to do with `2to3`).
@endtodo

Running the blocking analysis as
```bash
path/to/neci/utils/blocking.py -p 'plots/blocking.png' -f <numiter> -d24 -d23 -o/ > stats
```
will output a blocking plot to the plots subdirectory, starting after `<numiter>` iterations, which should be chosen at a point where the plateaus in plots (5) and (6) (i.e. the numerator and denominator for the error estimate) are both stable. In this example, we might choose `<numiter>` to be 9000. Running this, we get the following plot.

![](|media|/plots1/blocking.png)

Consisting of only three points, and having no plateau, this indicates that *we have not yet converged our FCIQMC calculation reliably.*
That is, if all the above 6 plots indicate convergence but the blocking analysis has no plateau (as in this example), it is most likely that you must continue the calculation to get more data.

### Continuing a NECI Calculation

In order to continue a NECI calculation (for example, if like in this example you have done a calculation only to find you do not have enough data), simply take the same NECI input as above, but add into the calc the `readPops` command, which indicates that NECI must read the popsfile previously created. You may also wish to increase the number of iterations `nmcyc`, e.g.
```text
title
  ...
  calc
  readPops
  nmcyc 70000
  ...
end
```
This will add data into the previous FCIMCStats. After you have run this the same way as described above, repeat the blocking analysis. The plots will all still have well-defined plateaus, but the blocking analysis will result in something like this:

![](|media|/plots2/blocking.png)

This time we see a clear plateau. The last point in this example indicates an excessively large blocking length (resulting in the estimate of the standard having too much noise). This is expected behaviour (but only in the blocking analysis; not in the other plots); what is important is that we see a plateau before the onset of this noise. Now, we may be more confident in our FCIQMC calculation. Above, we captured the output of the blocking script into a file called `stats`. In this example, the contents of that file is below:
```text
# of blocks mean (X_24)    std.err. (X_24) std.err.err. (X_24)  mean (X_23)    std.err. (X_23) std.err.err. (X_23)  cov(X_23,X_24) mean (X_24/X_23) std.err. (X_24/X_23)
722         -670.049359811 1.99777332e-01  5.26094278e-03       37244.1867452  9.82441101e+00  2.58716360e-01       -1.32220e+03   -0.017990709917  1.946944746518e-06
361         -670.049359811 2.68876787e-01  1.00204462e-02       37244.1867452  1.36690165e+01  5.09414168e-01       -1.28123e+03   -0.017990709917  1.911297631965e-06
180         -670.055219406 3.59855874e-01  1.90189739e-02       37244.4116806  1.86756159e+01  9.87036972e-01       -1.18406e+03   -0.017990758591  2.026001990428e-06
90          -670.055219406 4.67644101e-01  3.50514073e-02       37244.4116806  2.45252068e+01  1.83824197e+00       -1.01830e+03   -0.017990758591  2.124344681041e-06
45          -670.055219406 5.83345340e-01  6.21848221e-02       37244.4116806  3.08250668e+01  3.28596316e+00       -8.00932e+02   -0.017990758591  2.312712518500e-06
22          -670.085373646 7.34544351e-01  1.13342654e-01       37245.7631108  3.89686474e+01  6.01299283e+00       -6.25139e+02   -0.017990915414  2.494555852939e-06
11          -670.085373646 7.12029706e-01  1.59214682e-01       37245.7631108  3.79742226e+01  8.49129431e+00       -2.96473e+02   -0.017990915414  1.687704717679e-06
5           -669.988620008 6.61525494e-01  2.33884581e-01       37240.7990625  3.68826395e+01  1.30399822e+01       -1.21293e+02   -0.017990715475  1.907691051811e-06
2           -670.428455433 7.74026533e-02  5.47319410e-02       37267.1565039  9.84779297e+00  6.96344119e+00       -1.52449e+00   -0.017989793650  2.676810345979e-06
```

We wish to take the energy from the **first** row here in the second-to-last column and its corresponding uncertainty in the last column, i.e.:
```text
-0.017990709917 +/- 1.946944746518e-06
```

This is then our estimate for the correlation energy. To get the *total* energy, we must also add the reference energy, which can be found in the standard output of NECI (we called it `n2_neci.out`):
```text
Reference Energy set to:      -108.9606713172
```
(search for "Reference Energy"). You'll also find estimates for the correlation energy in the output file. However, this is not as trustworthy as doing a full blocking analysis.

### Final Steps

To be completely sure of our FCIQMC calculation, we must again continue from the popsfile with the `readPops` keyword, but change the number of walkers. The goal here is to verify that we have a sufficient number of walkers, and that we are converged with respect to the total number of walkers. Thus, to be really sure of our energy calculation, we must repeat the FCIQMC calculation but varying the number of walkers. The easiest way to do this is to restart from the previous popsfile and increase the total number of walkers. However, since the previous total number of walkers has already been reached, NECI is in variable-shift (or constant-walker-number) mode, and hence we need to tell NECI vary the number of walkers again.

To do this, we keep the `readPops` keyword and add the keyword and add the `walkContGrow` keyword into the `calc` block. We will, of course, also want to increase the total number of walkers (say, to 100000), and from our previous experience above we know we need more data so we can also increase the number of iterations, i.e.
```text
title
...
calc
# continue on from previous calculation
readPops
# allow growth from previous calc
walkContGrow

nmcyc 100000
...
end
```

Repeat this as above, do the same convergence analysis as above. Note, however, that since the number of iterations from which to start the blocking analysis (`<numiter>`) will be higher, as you will see by checking the plots (this is just because NECI needs some time to increase the walker number to the new target number and then stabilise).

Once you have done that, you may be much more confident about your calculated correlation energy.

Congratulations on your first FCIQMC calculation with NECI. The software has many more sophisticated options and can be used for bigger systems. The rest of this documentation discusses these in some details, though not in a tutorial format.
