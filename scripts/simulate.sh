#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00,h_data=8G
#$ -o ./job_out

PLINK_DIR=/u/project/pasaniuc/pasaniucdata/admixture/projects/admix-prs-uncertainty/data/PLINK
cur_dir="$(pwd)"

#chr_i=$((SGE_TASK_ID))
chr_i=$1
h2=$2
cau_prop=$3

mkdir -p /out/chr{$chr_i}_h{$h2}_cau{$cau_prop}/

python3 simulate.py \
    --data_dir ${PLINK_DIR} \
    --chr_i ${chr_i} \
    --h2 ${h2} \
    --cau_prop ${cau_prop} \
    --out ./out/chr${chr_i}_h2_${h2}_cau_${cau_prop}/

