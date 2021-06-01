# gnuplot script for visualising potential energy curves


# variables
m=0
nmax=5
lmax=6
elvl=2

# file settings
set datafile separator ","

filename(m,p,l,n) = sprintf(\
  "extracted_data/pot/m-%i.parity-%i.l_max-%i.n_basis_l_const-%i.txt",\
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

# parity = +1
parity=1

# vary N_{l}
set size 0.5, 0.5
set origin 0.0, 0.5
set title sprintf("Varying N_{l} [parity = %i, l_{max} = %i]", parity, lmax)
plot for [i=0:nmax] filename(m,parity,lmax,i) u 1:elvl

# vary l_{max}
set size 0.5, 0.5
set origin 0.5, 0.5
set title sprintf("Varying l_{max} [parity = %i, N_{l} = %i]", parity, 2**nmax)
plot for [i=0:lmax] filename(m,parity,i,nmax) u 1:elvl

# parity = -1
parity=-1

# vary N_{l}
set size 0.5, 0.5
set origin 0.0, 0.0
set title sprintf("Varying N_{l} [parity = %i, l_{max} = %i]", parity, lmax)
plot for [i=0:nmax] filename(m,parity,lmax,i) u 1:elvl

# vary l_{max}
set size 0.5, 0.5
set origin 0.5, 0.0
set title sprintf("Varying l_{max} [parity = %i, N_{l} = %i]", parity, 2**nmax)
plot for [i=0:lmax] filename(m,parity,i,nmax) u 1:elvl

pause -1 "enter any key to exit"
