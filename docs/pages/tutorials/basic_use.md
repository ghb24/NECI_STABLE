title: NECI Basic Use
author: Philip Haupt
---

[TOC]

# Basic NECI Tutorial

FCIQMC is not a blackbox method and as such may be daunting to first approach. In addition, there are many variants which excel in different kinds of problems. NECI has tools to make working with FCIQMC an easier task. While NECI is written predominantly in Fortran, in order to use it you need not know any programming language. In this tutorial, we assume an at least elementary understanding on the FCIQMC algorithm.[^fciqmc]

The goal of this tutorial is to provide a practical, brief supplement to the full [NECI user's guide](../01_neci_user.html). We will use NECI to calculate the ground-state energy of the Nitrogen dimer in a (6,6) active space. We use the equilibrium geometry, 2.074 bohr radii.

@note
For this problem, FCIQMC will actually underperform compared to conventional Davidson CI. In addition, Davidson CI has no stochastic error (unlike FCIQMC), so it is strongly preferred for this problem. However, this is a relatively cheap example with which to get comfortable with NECI and FCIQMC. FCIQMC has better scaling behaviour, can be better parallelised, and benefits from sparsity. If you go to (e.g.) an (18,18) CAS then FCIQMC outshines Davidson, which is not a feasible method for problems of this size.
@endnote

First, we must generate the FCIDUMP file. NECI is *not* a standalone quantum chemistry suite, and cannot do this itself. For this, choose any program that can generate these (e.g. PySCF, Molpro, Molcas). This has been done for you, and you may download the file [here](FCIDUMP).

@todo
Perhaps also include the Molpro code to generate this FCIDUMP file?
@endtodo

## Anatomy of an Input File

In order to run a NECI calculation, we must create an input file (extension `.inp`). Here is an example for the FCIDUMP file provided, called `n2_neci.inp`.
```
# comments are given like this
( or like this

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

All these keywords (and plenty more) are explained in the [NECI user's guide](../01_neci_user.html). Let's break down the structure of this input file.

- Comments can be added in the code with either `#` or `(` at the start of a line.
@todo
I was told that only `(` is supposed to work but I have been using `#`. May this cause some kind of problems later, or may I keep `#`?
@endtodo
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
    - Another form of FCIQMC is i-FCIQMC, which uses "initiator" walkers to speed up the calculation.[^initiator] This is optional, but generally recommended. Presence of the `truncinitiator` keyword in the `calc` block indicates i-FCIQMC, and `addToInitiator 3` means that any determinant with a population >= 3 will become an initiator.
    - Finally, as a mandatory subblock *inside* the calc block, we have the `methods` subblock, which necessarily ends with `endmethods`. Here, we simply specify what kind of calculation to do. We choose `method vertex fcimc`, which simply means to run an FCIQMC calculation.
- Finally, we have an *optional* `logging` block which ends with `endlog`. By default, NECI will keep track of the population of walkers in a "POPSFILE", which is by default ASCII. Here, we wish to have an HDF5 POPSFILE which is generally better-performing. To do this, inside the `logging` block we have the `hdf5-pops` keyword.

@note
You cannot use the `hdf5-pops` keyword if you did not build NECI with HDF5.
@endnote

## Running NECI

After building, the NECI executable will be in `path/to/neci/build/bin/neci` (e.g. if I installed NECI in my home directory, it would be `~/neci/build/bin/neci`).

NECI must be run with MPI (there is no serial version of NECI). To run the above input file, we must run
```bash
mpirun -np 4 path/to/neci/build/bin/neci n2_neci.inp
```
where you can replace the `4` with however many processors with which you want to run (4 being a *very* modest number). This will print a lot of standard output, which you may wish to capture, e.g.
```bash
mpirun -np 4 ~/neci/build/bin/neci n2_neci.inp > n2_neci.out
```

@note
If you made a mistake somewhere in your file (for example, a typo), NECI tries to give useful error messages. However, MPI can sometimes interfere with this. Before doing a full calculation, you may want to run simply without `mpirun`, e.g. simply
```
~/neci/build/bin/neci n2_neci.inp > n2_neci.out
```
If this starts to run without an error, then you may stop it and run with `mpirun`.
@endnote

@note
NECI has a useful utility to dynamically change variables *while it is running*. To do this, create a file called `CHANGEVARS`, and input whatever keyword you wish to change, e.g. if you think `nmcyc 10000` is too small, enter into `CHANGEVARS` the text `nmcyc 20000` (or whatever else). Then, the value will be updated in the next iteration of NECI. This is helpful for interacting with a running simulation.
@endnote

@warning
It is generally not advisable to terminate a program with CTRL+C. Instead, use the above `CHANGEVARS` trick, to create a soft exit. Simple open/create `CHANGEVARS`, type in `SOFTEXIT` and save. Alternatively, you can do this in one line on the terminal with
```
echo SOFTEXIT > CHANGEVARS
```
@endwarning

This calculation should produce a few other files in the directory, namely:
- `popsfile.h5` (or `POPSFILE` if you did not include `hdf5-pops`)
- `INITATORStats`
- `FCIMCStats`

The popsfile will be very useful in case we wish to continue running the simulation. `FCIMCStats` has columns of useful data, which we will explore now.

## Checking Convergence and Analysing Results

The file `FCIMCStats` has several useful columns which you will want to plot to ensure convergence. To do this in one line, [there is a convenience script here](plot_fcimcstats.plt), which is run with gnuplot via `gnuplot plot_fcimcstats.plt` and will output plots of the most useful columns to a new `plots/` directory. Your results should look something like this:

@todo
these plots are not showing up -- maybe need media_dir keyword
@endtodo

![](|media|/plots1/check_totE.png)

<!-- ![](plots1/check_totE.png)
![](plots1/check_totWalkers.png)
![](plots1/check_denominator.png)
![](plots1/check_numerator.png)
![](plots1/check_refPop.png)
![](plots1/check_shift_energy.png) -->



@todo
![](plots1/blocking.png)
@endtodo

@todo
mention dynamically changing stuff and soft exit with CHANGEVARS.
mention checking input by running _without_ mpirun.
mention restarts
spend time on convergence analysis
@endtodo

[^fciqmc]: Booth, G. H., Thom, A. J. W. & Alavi, A. Fermion monte carlo without fixed nodes: A game of life, death, and annihilation in Slater determinant space. J. Chem. Phys. 131, 054106 (2009).
[^initiator]: D. Cleland, G.H. Booth, A. Alavi, J. Chem. Phys. 132, 041103 (2010)