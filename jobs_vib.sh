#!/bin/bash

# constant parameters
alpha_l_const="1.0"
d_r="0.01"
r_max="10.0"

# parameter sets
reduced_mass_set="\
918.07635 1223.89925 1376.39236 1835.24151 2200.87999 2748.46079"
n_basis_l_const_set="1 2 4 8 16 32 64 128"

# compile
make vibrational

# run jobs
echo "H2+ vibrational wave function jobs"
echo "parameter constants:"
echo "> alpha_l_const: ${alpha_l_const}"
echo "> d_r: ${d_r}"
echo "> r_max: ${r_max}"
echo "parameter sets: "
echo "> reduced_mass_set: ${reduced_mass_set}"
echo "> n_basis_l_const_set: ${n_basis_l_const_set}"
echo "> "

echo "execution:"
initial=$(date +%s)
job=1
for n_basis_l_const in ${n_basis_l_const_set} ; do
  for reduced_mass in ${reduced_mass_set} ; do
    start=$(date +%s)
    printf "%4i @ [%is], " ${job} $((start - initial))

    printf "params (%3i, %10.5f), " \
           ${n_basis_l_const} ${reduced_mass}

    bin/vibrational  ${n_basis_l_const} ${alpha_l_const} ${reduced_mass} \
                     ${d_r} ${r_max}

    end=$(date +%s)
    runtime=$((end - start))
    printf "runtime [%is]\n" ${runtime}

    ((job=job+1))
  done
done
