#!/bin/bash

eval $(spack load --sh singularityce@3.8.0)

img=$1

run_cmd=$2

singularity exec \
  -B /scratch/mblab/chasem \
  -B /ref/mblab/data \
  -B "$PWD" \
  $img \
  /bin/bash -c "cd $PWD; $run_cmd"
