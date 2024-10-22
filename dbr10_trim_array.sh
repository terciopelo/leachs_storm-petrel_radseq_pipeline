#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 4-0:0:0
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --job-name=ParseFastQ_DBR10_batch2array
#SBATCH --array=10-17
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020 python/2.7.18

python ./ParseDBR_ddRAD/ParseFastQ.py -r x${SLURM_ARRAY_TASK_ID}_R1.fastq  \
         -R x${SLURM_ARRAY_TASK_ID}_R2.fastq \
         -i TAGCTT -e AATT -n ./batch2/batch2_DBR10_x${SLURM_ARRAY_TASK_ID}_R1_output.fastq \
         -N ./batch2/batch2_DBR10_x${SLURM_ARRAY_TASK_ID}_R2_output.fastq \
         --drop ./batch2/batch2_dropped_DBR10_${SLURM_ARRAY_TASK_ID}.txt -Z -l 2
