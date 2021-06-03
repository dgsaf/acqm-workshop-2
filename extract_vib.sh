#!/bin/bash

# extraction parameters
f_max="20"
scale="0.0035"

# constant parameters
alpha_l_const="1.0"
d_r="0.01"
r_max="100.0"

# directory functions
# - ${1} = ${n_basis_l_const}
# - ${2} = ${reduced_mass}
job_dir () {
  local fmt="\
output/vib/\
n_basis_l_const-%i.\
alpha_l_const-%.5f.\
reduced_mass-%.5f.\
d_r-%.5f.\
r_max-%.5f"
  printf "${fmt}" ${1} ${alpha_l_const} ${2} ${d_r} ${r_max}
}

extract_file () {
  local fmt="\
extracted_data/vib/\
n_basis_l_const-%i.\
reduced_mass-%.5f.shifted_radial.txt"
  printf "${fmt}" ${1} ${2}
}

# parameter sets
reduced_mass_set="\
  918.07635 1223.89925 1376.39236 1835.24151 2200.87999 2748.46079"
n_basis_l_const_set="1 2 4 8 16 32 64 128"

# run jobs
echo "H2+ vibrational wave function extraction"
echo "parameter constants:"
echo "> alpha_l_const: ${alpha_l_const}"
echo "> d_r: ${d_r}"
echo "> r_max: ${r_max}"
echo "parameter sets: "
echo "> reduced_mass_set: ${reduced_mass_set}"
echo "> n_basis_l_const_set: ${n_basis_l_const_set}"
echo "> "

echo "extraction:"
job=1
for n_basis_l_const in ${n_basis_l_const_set} ; do
  for reduced_mass in ${reduced_mass_set} ; do
    # determine output directory for these parameters
    outdir=$(job_dir ${n_basis_l_const} ${reduced_mass})

    # ensure output directory exists
    if [ ! -d ${outdir} ] ; then
      printf "%i @ %s does not exist\n" ${job} ${outdir}
    else
      printf "%i @ extracting wavefunctions from %s\n" ${job} ${outdir}

      energyfile="${outdir}/eigen_values.txt"
      # printf "%i @ energy file %s\n" ${job} ${energyfile}

      radialfile="${outdir}/eigen_radial.txt"
      # printf "%i @ radial file %s\n" ${job} ${radialfile}

      extfile=$(extract_file ${n_basis_l_const} ${reduced_mass})
      # printf "%i @ extract file %s\n" ${job} ${extfile}

      # extract energies
      energies="0.0"
      for (( k = 1 ; k<=${n_basis_l_const} && k<=${f_max} ; k++ )) ; do
        energy=$(sed -n "${k}p" ${energyfile})
        energies="${energies} ${energy}"

        printf "%i , %i : %f\n" ${n_basis_l_const} ${k} ${energy}
      done
      echo "${energies}"

      # shift radial wave functions by energy
      tempfile="${extfile}.temp"
      echo "${energies}" > ${tempfile}
      cat ${radialfile} >> ${tempfile}

      awk -v FM="${f_max}" -v SC="${scale}" \
          ' {\
              if (FNR==1) {\
                for (i=1; i<=NF && i<=FM+1; i++) {\
                  en[i]=$i;\
                }\
              }\
              else {\
                printf "%f ", $1;\
                for (i=2; i<=NF && i<=FM+1; i++) {\
                  shifted = (SC*$i) + en[i];\
                  printf "%f ", shifted;\
                }\
                printf "\n";\
              }\
            }'\
          ${tempfile} > ${extfile}

      rm ${tempfile}
    fi
    ((job=job+1))
  done
done
