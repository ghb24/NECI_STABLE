title: Tutorials 
---

TODO : come up with tutorials, including example input and going over what they mean. Also (most important imo) go over the convergence checks and error analysis, as applied to a real problem (e.g. N2? Maybe even just STO-3G H2). Include a case where you need to restart a calculation because NECI wasn't fully converged. Go over the output file as well 

TODO some kind of documentation for the columns of `FCIMCStats`.

TODO discuss what specifically should be put here 

also here is a simple gnuplot script to do all the relevant plots (after this you still need to call `blocking.py`): 

```
# plot_fcimcstats.plt
system("mkdir plots")

datafile = "FCIMCStats"
outdir = "plots/"
set terminal png # size 500,500

# want to plot cols 5, 11(?), 12, 24, 25, maybe 2?
# maybe some output instructions? 

# col 5 Total Walkers 
outfile = outdir."check_totWalkers.png" 
set output outfile 
set xlabel "Step" 
set title "Total Walkers"
p "FCIMCStats" u 1:5 w lp 

# col 12 Reference population 
outfile = outdir."check_refPop.png" 
set output outfile 
set xlabel "Step" 
set title "Population on Reference Determinant"
p "FCIMCStats" u 1:12 w lp 

# col 24 denominator
outfile = outdir."check_denominator.png" 
set output outfile 
set xlabel "Step" 
set title "Denominator (Col 24)"
p "FCIMCStats" u 1:24 w lp 

# col 25 numerator 
outfile = outdir."check_numerator.png" 
set output outfile 
set xlabel "Step" 
set title "Numerator (Col 25)"
p "FCIMCStats" u 1:25 w lp 

# Comparison plot 
outfile = outdir."check_shift_energy.png"
set output outfile
# set xrange [400000:1000000]
set yrange [:*<0]
p "FCIMCStats" u 1:2 w lp # , "FCIMCStats" u 1:9 w lp
MAXY=GPVAL_Y_MAX
MINY=GPVAL_Y_MIN
MAXX=GPVAL_X_MAX
MINX=GPVAL_X_MIN
unset autoscale
set yrange [MINY:MAXY]
set xrange [MINX:MAXX]
set output outfile
set xlabel "Step" 
set ylabel "Energy Estimate"
set title "compare col 11,2,25/24"
p "FCIMCStats" u 1:($25/$24) w lp t "col 25/col 24", "FCIMCStats" u 1:11 w lp t "Shift", "FCIMCStats" u 1:2 w lp t "Proj.E"
```