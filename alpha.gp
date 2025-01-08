set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'alpha_rich.png'

set title "Convergence of Richardson (optimal alpha)"
set xlabel "Iteration"
set ylabel "Residual Norm"

set grid
#set logscale y

plot "ALPHA_RESVEC.dat" using 0:1 with linespoints title "Richardson alpha"