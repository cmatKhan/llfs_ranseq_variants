#!/usr/bin/bash

#SBATCH --mem-per-cpu=100M
#SBATCH -J pull_from_s3
#SBATCH -o pull_from_s3.out

# usage:
# 
# sbatch \
#   --array=1-194 
#   path/to/scripts/pull_from_s3.sh \
#   lookups/20230617.txt \
#   data/20230617
#

eval $(spack load --sh py-s3cmd@2.3.0)

read s3Path data_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p $1)

mkdir -p ${data_dir}

s3cmd get ${s3Path} ${data_dir}/
