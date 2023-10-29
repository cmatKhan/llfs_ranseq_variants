#!/bin/bash

#SBATCH --mem-per-cpu=2G
#SBATCH -N 1
#SBATCH -n 11
#SBATCH --time=02:00:00
#SBATCH -J compile_results
#SBATCH -o compiles_results_%A.out

eval $(spack load --sh singularityce@3.11.4)

singularity_image=/scratch/mblab/chasem/software/llfs_rnaseq_latest.sif

dir=$1

singularity exec \
  -B /scratch/mblab \
  -B /ref/mblab \
  -B "$PWD" \
  $singularity_image \
  /bin/bash -c \
  "cd $PWD; \
   export R_LIBS=/project/renv/library/R-4.2/x86_64-pc-linux-gnu; \
   /scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_variants/scripts/compile_compare_results.R --dir "$dir" "

