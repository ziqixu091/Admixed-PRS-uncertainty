#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00,h_data=8G
#$ -o ./job_out

PLINK_DIR=/u/project/pasaniuc/pasaniucdata/admixture/projects/admix-prs-uncertainty/data/PLINK

#chr_i=$((SGE_TASK_ID))
chr_i=$1
var_g=$2
var_e=$3

python3 simulate.py \
    --data_dir ${PLINK_DIR} \
    --chr_i ${chr_i} \
    --var_g ${var_g} \
    --var_e ${var_e} 

