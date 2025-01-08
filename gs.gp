set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'gauss_seidel.png'

set title "Convergence of Gauss-Seidel Method"
set xlabel "Iteration"
set ylabel "Residual Norm"
set grid
#set logscale y

plot "GS_RESVEC.dat" using 0:1 with linespoints title "Gauss-Seidel"