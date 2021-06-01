# gnuplot script for visualising potential energy curves

# file settings
set datafile separator ","

filename(m,p,l,n) = sprintf(\
  "extracted_data/m-%i.parity-%i.l_max-%i.n_basis_l_const-%i.txt",\
  m, p, l, 2**n)

# common settings
set xlabel "rz"
set ylabel "energy"
unset key
set style data linespoints

# multiplot 2x2
set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot

# parity = +1, vary N_{l}
set size 0.5, 0.5
set origin 0.0, 0.5
set title "Varying N_{l} [parity = +1, l_{max} = 4]"
plot for [i=0:5] filename(0,1,4,i) w lines

# parity = +1, vary l_{max}
set size 0.5, 0.5
set origin 0.5, 0.5
set title "Varying l_{max} [parity = +1, N_{l} = 32]"
plot for [i=0:4] filename(0,1,i,5) w lines

# parity = -1, vary N_{l}
set size 0.5, 0.5
set origin 0.0, 0.0
set title "Varying N_{l} [parity = -1, l_{max} = 4]"
plot for [i=0:5] filename(0,-1,4,i) w lines

# parity = -1, vary l_{max}
set size 0.5, 0.5
set origin 0.5, 0.0
set title "Varying l_{max} [parity = -1, N_{l} = 32]"
plot for [i=0:4] filename(0,-1,i,5) w lines

pause -1 "enter any key to exit"
