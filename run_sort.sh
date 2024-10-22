#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=sort-lesp-bams
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --mem 48G
#SBATCH -c 1
#SBATCH --array=10
#SBATCH --time 4:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load samtools

for FILE in `cat x${SLURM_ARRAY_TASK_ID}`; do
	samtools sort -o ./stacks_align/${FILE}_sorted.bam ./stacks_align/${FILE}.bam
	samtools index ./stacks_align/${FILE}_sorted.bam
done
