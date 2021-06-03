# gnuplot script for visualising H2+ vibrational wave functions


# visualisation constants
nf_max=10

# constants
alpha=1.0
d_r=0.05
r_max=10.0
n2_max=8

# file settings
potential = sprintf("analytic_data/PEC.1ssg")

directory(n, mu) = sprintf(\
  "output/vib/n_basis_l_const-%i.alpha_l_const-%.5f.reduced_mass-%.5f.d_r-%.5f.r_max-%.5f",\
  n, alpha, mu, d_r, r_max)

extracted(n, mu) = sprintf(\
  "extracted_data/vib/n_basis_l_const-%i.reduced_mass-%.5f.shifted_radial.txt",\
  n, mu)

# common settings
set xlabel "rz"
set ylabel "energy"
unset key
set style data lines

set xrange [0:6]
set yrange [-1.5:1.5]

# H2
mu=918.07635

do for [n2=1:n2_max] {
  n = 2**n2

  plot potential

  n_max = (n < nf_max) ? n : nf_max
  do for [i=2:n_max] {
    plot potential, for [j=2:i] extracted(n, mu) u 1:j
    pause 0.5
  }
  pause 1.0
}
pause -1 "enter any key to exit"
