#!/bin/bash

# usage 
# sbatch --array=1-1408 --mem=3G --cpus-per-task=1 --time=20 compare_samples.sh sex_mislabel_compare_redo.tsv

# $1 is the lookup

#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=1
#SBATCH -J rna_dna_compare
#SBATCH -o rna_dna_compare.out

eval $(spack load --sh singularityce@3.11.4)

singularity_image=/scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_latest.sif

read rna_subject_id visit dna_subject_id rna_gds chr dna_gds < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")

echo "$chr"

singularity exec \
  -B /scratch/mblab \
  -B "$PWD" \
  -B /ref/mblab/data \
  $singularity_image \
  /bin/bash -c \
  "cd $PWD; \
   /scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_variants/scripts/extract_rna_dna_overlap_genotypes.R \
   --chr $chr \
   --rna_sample $rna_subject_id \
   --rna_visit $visit \
   --dna_sample $dna_subject_id \
   --dna $dna_gds \
   --rna $rna_gds"
