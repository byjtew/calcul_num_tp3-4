set grid

set ylabel "Nb of iterations to converge"

set xlabel "Matrix leading dimension"

plot "jacobi_gauss_seidel_richardson_results.dat" using 1:2 with linespoints title "Jacobi", "jacobi_gauss_seidel_richardson_results.dat" using 1:3 with linespoints title "Gauss-Seidel", "jacobi_gauss_seidel_richardson_results.dat" using 1:4 with linespoints title "Richardson"