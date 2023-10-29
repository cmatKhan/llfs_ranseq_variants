#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=20
#SBATCH -J rna_dna_comp.out
#SBATCH -o rna_dna_comp.out

# $1 is the lookup

# usage: sbatch --array=1-98 ../scripts/conduct_s3_compare.sh ../lookups/chunk_9_lookup.txt

# exit immediately if error. Note that this is 
# considered bad practice. Preferred method is 
# trap 'do something' ERR
set -e

eval $(spack load --sh singularityce@3.11.4)
eval $(spack load --sh py-s3cmd@2.3.0)

singularity_image=/scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_latest.sif

read s3_vcf s3_vcf_index rna_id subject_id visit chr dna_gds < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")

echo "pulling $(basename $s3_vcf) from s3"
s3cmd get --skip-existing $s3_vcf .
s3cmd get --skip-existing $s3_vcf_index .
echo "done pulling $(basename $s3_vcf) from s3"

echo "converting $(basename $s3_vcf) to gds"
vcf=$(basename $s3_vcf)
vcf_index=$(basename $s3_vcf_index)


singularity exec \
  -B /scratch/mblab \
  -B /ref/mblab \
  -B "$PWD" \
  $singularity_image \
  /bin/bash -c \
  "cd $PWD; \
   export R_LIBS=/project/renv/library/R-4.2/x86_64-pc-linux-gnu; \
   /scratch/mblab/chasem/llfs_rna_dna_compare_test/llfs_rnaseq_variants/R/convert_vcf_to_gds.R \
   --vcf $vcf"
echo "done with converting $(basename $s3_vcf) to vcf" 

echo "sending $(basename $s3_vcf) to s3" 
s3path_dirname=$(dirname $s3_vcf)
gds=$(basename ${s3_vcf%.vcf.gz}.gds)

s3cmd sync $gds $s3path_dirname/
echo "done sending $(basename $s3_vcf) to s3"

echo "comparing rna against dna: $(basename $s3_vcf)"
visit=$visit
rna_sample=$rna_id
rna=$gds
dna_sample=$subject_id
dna=$dna_gds

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
   --dna $dna \
   --rna $rna"

rm $vcf
rm $vcf_index
rm $gds

echo complete
