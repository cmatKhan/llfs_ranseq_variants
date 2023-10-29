#!/bin/bash

#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-4556
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:60:00
#SBATCH -J rna_dna_comp
#SBATCH -o rna_dna_comp_%A.out

# usage:
# sbatch while read rna_id visit rna_gds output_dir; do ../scripts/only_compare.sh $rna_id $visit $rna_gds $dna_gds; done < ../samplesheets/test.txt
# where a line of the samplesheet looks like:
# 43319779	3	43319779_v3_T1.haplotypecaller.filtered.gds	/path/to/output/dir
# and has columns rna_gds_identifier (the identifier for that library in the gds file), visit, rna_gds, output_dir

eval $(spack load --sh singularityce@3.11.4)

singularity_image=/scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_latest.sif

wgs_dna_subjects=/scratch/mblab/chasem/llfs_rna_dna_compare_test/lookups/wgs_dna_subject_ids.txt
dna_gds=/ref/mblab/data/llfs/agds/LLFS.WGS.freeze5.chr1.gds
chr=1

rna_sample=$1
visit=$2
rna_gds=$3
output_dir=$4

read dna_sample < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$wgs_dna_subjects")

singularity exec \
  -B /scratch/mblab \
  -B /ref/mblab \
  -B "$PWD" \
  $singularity_image \
  /bin/bash -c \
  "cd $PWD; \
   export R_LIBS=/project/renv/library/R-4.2/x86_64-pc-linux-gnu; \
   /scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_variants/scripts/extract_rna_dna_overlap_genotypes.R \
   --chr $chr \
   --rna_sample $rna_sample \
   --rna_visit $visit \
   --dna_sample $dna_sample \
   --dna $dna_gds \
   --rna $rna_gds \
   --output_prefix $output_dir"

