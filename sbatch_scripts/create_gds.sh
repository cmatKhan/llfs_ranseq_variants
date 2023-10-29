#!/bin/bash

# $1 is the lookup

#SBATCH -J vcf_to_gds
#SBATCH -o vcf_to_gds.out

export SINGULARITY_TMPDIR=/tmp

eval $(spack load --sh singularityce@3.8.0)

IMAGE=/scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_latest.sif

read vcf < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")

singularity exec \
  -B /scratch/mblab \
  -B /ref/mblab/data \
  -B "$PWD" \
  $IMAGE \
  /bin/bash -c \
  "cd $PWD; \
   export R_LIBS=/project/renv/library/R-4.2/x86_64-pc-linux-gnu; \
   /scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_variants/R/convert_vcf_to_gds.R \
   --vcf $vcf"
