system("mkdir plots")

datafile = "FCIMCStats"
outdir = "plots/"
set terminal png # size 500,500

# col 5 Total Walkers
outfile = outdir."check_totWalkers.png"
set output outfile
set xlabel "Step"
set title "Total Walkers (Col 5)"
p "FCIMCStats" u 1:5 w lp

# col 12 Reference population
outfile = outdir."check_refPop.png"
set output outfile
set xlabel "Step"
set title "Population on Reference Determinant (Col 12)"
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

# col 23 total energy
outfile = outdir."check_totE.png"
set output outfile
set xlabel "Step"
set title "Total Energy (Col 23)"
p "FCIMCStats" u 1:23 w lp

# Comparison plot
outfile = outdir."check_shift_energy.png"
set output outfile
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
set title "Shift and Projected Energies"
p "FCIMCStats" u 1:2 w lp t "Shift", "FCIMCStats" u 1:11 w lp t "Proj.E"