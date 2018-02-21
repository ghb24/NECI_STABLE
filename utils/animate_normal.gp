reset
set term gif animate 10 size 1024,1024
set output "normal.gif"
set yrange[-0.1:0.1]
set xlabel "Det"
set ylabel "Amplitude"
#Line style of the instantaneous wavefunction
set style line 1 lt 4 lw 1 pt 2 ps 0.2
#Line style of the exact wavefunction
set style line 2 lt 1 lw 2
#Line style of the average wavefunction
set style line 3 lt 2 lw 2 pt 3 ps 0.2
set print "-"
do for [iter=1000:100000:1000]{
    print "Iteration: ".iter
    set title "Normal Initiator Method - Iteration: ".iter
    plot \
         'SpawnHist-'.iter u 1:4 w lp ls 1 t "Instant.",\
         'SymDETS'         u 1:4 w l  ls 2 t "Exact",\
         'SpawnHist-'.iter u 1:2 w lp ls 3 t "Avg."  
}
