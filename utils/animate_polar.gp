reset
set term gif animate 10 size 1024,1024
set output "polar.gif"

set polar
unset raxis
#unset key
unset rtics
unset border
unset xtics
unset ytics

#Miimum and Maximum Exciations
min_ex = 2
max_ex = 8 

#Exponential (as a function of excitation level) scaling factor for the population
scale = 2

#Distance between exciations
dist = 1

#Line style of the instantaneous wavefunction
set style line 1 lt 4 lw 1 pt 2 ps 0.2
#Line style of the exact wavefunction
set style line 2 lt 1 lw 2
#Line style of the average wavefunction
set style line 3 lt 2 lw 2 pt 3 ps 0.2

set xrange [-max_ex*dist:+max_ex*dist]
set yrange [-max_ex*dist:+max_ex*dist]

#Find the line index of the first det of an excaitation
index(ex) = system("awk '{ if (\$2 == ".ex.") {print \$1; exit }}' SymDETS")

set print "-"
do for [iter=100:1000:100]{
    print "Iteration: ".iter
    set multiplot
    set title "Normal Initiator Method - Iteration: ".iter
    #We need the exact solution, so we append the columns of SymDETS
    system("paste -d' ' SpawnHist-".iter." SymDETS > tmp")
    do for [ex=min_ex:max_ex]{
        #Find the first and last line for excitation level
        begin = index(ex)
        end = index(ex+1)-1

        #Extract the lines of an exciation
        system("sed -n '".begin.",".end."p' tmp > ex_tmp")
        #Repeat the first line at the end to close the polar graph
        system("sed -n '1p' ex_tmp  >> ex_tmp")

        x2angle(x) = (2*pi*(x-begin)/(end-begin+1))
        y2radius(y) = (ex*dist+(scale**(ex)*y))

        plot \
             'ex_tmp' u (x2angle($1)):(y2radius($4))  w lp ls 1 t "Instant.",\
             'ex_tmp' u (x2angle($1)):(y2radius($13)) w l  ls 2 t "Exact",\
             'ex_tmp' u (x2angle($1)):(y2radius($2))  w l ls 3 t "Avg."
        #smooth bezier
    }
    unset multiplot
}
set output
