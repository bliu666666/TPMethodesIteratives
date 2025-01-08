set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'jacobi.png'

set title "Convergence of Jacobi Method"
set xlabel "Iteration"
set ylabel "Residual Norm"
set grid
#set logscale y

plot "JAC_RESVEC.dat" using 0:1 with linespoints title "Jacobi"