# gnuplot script for visualising H2+ vibrational wave functions

# constants
alpha=1.0
nmax=5
d_r=0.1
r_max=100.0

# file settings
potential = sprintf("analytic_data/PEC.1ssg")

directory(n, mu) = sprintf(\
  "output/vib/n_basis_l_const-%i.alpha_l_const-%.4f.reduced_mass-%.4f.d_r-%.4f.r_max-%.4f",\
  n, alpha, mu, d_r, r_max)

radial(n, mu) = sprintf("%s/eigen_radial.txt", directory(n, mu))

# common settings
set xlabel "rz"
set ylabel "energy"
unset key
set style data lines

set xrange [0:10]
set yrange [-1.5:1.5]

# H2
mu=918.0000
n=10

plot potential
do for [i=2:n] {
  plot potential, for [j=2:i] radial(n, 918.0) u 1:j
  pause 1.0
}

pause -1 "enter any key to exit"
