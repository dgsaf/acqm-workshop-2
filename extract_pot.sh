#!/bin/bash

# constant parameters
alpha_l_const="1.0"
nuclei_charge="1"
lambda_max="10"
d_r="0.1"
r_max="75.0"
d_rz="0.5"
rz_max="10.0"

# parameter sets
m_set="0"
parity_set="-1 +1"
l_max_set="0 1 2 3 4 5 6"
n_basis_l_const_set="1 2 4 8 16 32"

# common part of all output directory paths
common=$(printf "alpha_l_const-%.4f.nuclei_charge-%i.lambda_max-%i" \
                ${alpha_l_const} ${nuclei_charge} ${lambda_max})
common=$(printf "%s.d_r-%.4f.r_max-%.4f/" ${common} ${d_r} ${r_max})

# pattern to match axial distance sub-directories on
pattern="rz-([0-9].+)"

# extract jobs
job=1
for n_basis_l_const in ${n_basis_l_const_set} ; do
  for l_max in ${l_max_set} ; do
    for parity in ${parity_set} ; do
      for m in ${m_set} ; do
        specific=$(printf "output/pot/m-%i.parity-%i" \
                          ${m} ${parity})
        specific=$(printf "%s.l_max-%i.n_basis_l_const-%i" \
                          ${specific} ${l_max} ${n_basis_l_const})

        outdir="${specific}.${common}"

        if [ ! -d ${outdir} ] ; then
          printf "%i @ %s does not exist\n" ${job} ${outdir}
        else
          printf "%i @ extracting energies from %s\n" ${job} ${outdir}

          # construct extract file
          extractfile=$(printf "extracted_data/pot/m-%i.parity-%i" \
                               ${m} ${parity})
          extractfile=$(printf "%s.l_max-%i.n_basis_l_const-%i.txt" \
                               ${extractfile} ${l_max} ${n_basis_l_const})

          # pattern match for axial distance sub-directories
          subdirs=$(ls ${outdir})

          for subdir in ${subdirs} ; do
            if [[ $subdir =~ $pattern ]] ; then
              rz=${BASH_REMATCH[1]}

              printf "%i @ extracting energy for rz=%s\n" ${job} ${rz}

              energyfile="${outdir}rz-${rz}/eigen_values.txt"

              energy=$(awk -vORS=, '{ print $1 }' ${energyfile} | sed 's/,$//')
              # printf "%i @ %s, %s \n" ${job} ${rz} ${energy}

              printf "%s,%s\n" ${rz} ${energy} >> ${extractfile}
            else
              printf "%i @ %s is not an axial-distance sub-directory \n" \
                     ${job} ${subdir}
            fi
          done

          # if extract file exists, sort it on axial distance
          if [ -f ${extractfile} ] ; then
            printf "%i @ sorting extract file\n" ${job}

            sort -o ${extractfile} -k1 -n ${extractfile}
          fi

        fi

        printf "%i @ extracted\n" ${job}

        ((job=job+1))
      done
    done
  done
done

printf "extraction complete\n"
