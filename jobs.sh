#!/bin/bash

# constant parameters
alpha_l_const="1.0"
nuclei_charge="1"
lambda_max="20"
d_r="0.1"
r_max="100.0"
d_rz="0.25"
rz_max="5.0"

# parameter sets
m_set="0"
parity_set="-1 1"
l_max_set="0 1 2 3 4"
n_basis_l_const_set="1 2 4 8 16 32 64 128"

# compile
make potential_curves

# run jobs
echo "H2+ potential energy curve jobs"
echo "parameter constants:"
echo "> alpha_l_const: ${alpha_l_const}"
echo "> nuclei_charge: ${nuclei_charge}"
echo "> lambda_max: ${lambda_max}"
echo "> d_r: ${d_r}"
echo "> r_max: ${r_max}"
echo "> d_rz: ${d_rz}"
echo "> rz_max: ${rz_max}"
echo "parameter sets: "
echo "> m_set: ${m_set}"
echo "> parity_set: ${parity_set}"
echo "> l_max_set: ${l_max_set}"
echo "> n_basis_l_const_set: ${n_basis_l_const_set}"

for m in ${m_set} ; do
    for parity in ${parity_set} ; do
        for l_max in ${l_max_set} ; do
            for n_basis_l_const in ${n_basis_l_const_set} ; do
                bin/potential_curves ${m} ${parity} ${l_max} \
                                     ${n_basis_l_const} ${alpha_l_const} \
                                     ${nuclei_charge} ${lambda_max} \
                                     ${d_r} ${r_max} ${d_rz} ${rz_max}

            done
        done
    done
done
