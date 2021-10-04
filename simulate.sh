#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00,h_data=16G,highp -pe shared 1
#$ -o ./job_out

. /u/local/Modules/default/init/modules.sh
module load python/3.7.0

cur_dir = "$(pwd)"
mkdir -p sim/
PLINK_DIR = "/u/project/pasaniuc/pasaniucdata/admixture/projects/admix-prs-uncertainty/data/PLINK"

#chr_i=$((SGE_TASK_ID))
chr_i=$1
var_g=$2
var_e=$3

python3 simulate.py \
    --data_dir ${PLINK_DIR} \
    --chr_i ${chr_i} \
    --var_g ${var_g} \
    --var_e ${var_e} \
    --out_dir ${cur_dir}/sim

