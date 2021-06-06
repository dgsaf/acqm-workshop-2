# H2+ Molecular Structure

## Contents
  - `report/`
  - `src/`
  - `obj/`
  - `mod/`
  - `bin/`
  - `output/`
  - `extracted_data/`
  - `analytic_data/`
  - `makefile`
  - `jobs_pot.sh`
  - `extract_pot.sh`
  - `visualise_pot.sh`
  - `jobs_vib.sh`
  - `extract_vib.sh`
  - `visualise_vib.sh`

## Potential Curves
  1. Modify `makefile` to change compiler, optimisation level etc.
  1. Synchronise the job settings across `jobs_pot.sh`, `extract_pot.sh` and
     `visualise_pot.plt`.
  2. Run `bash jobs_pot.sh` to write H2+ data to `output/pot/` across a range of
     axial distance values. Note that if a calculation for a given set of
     parameters has already been performed, it will be skipped over.
  3. Run `bash extract_pot.sh` to extract the potential energy curve for the job
     set to `extracted_data/pot/`.
  4. Run `gnuplot visualise_pot.plt` for a simple visualisation.

## Vibrational
Same as for [the previous section](#user-content-potential-curves), except with
`*_vib` rather than `*_pot` for script files.
