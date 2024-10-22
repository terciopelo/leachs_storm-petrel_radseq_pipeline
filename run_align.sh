#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=align_lesp_rad
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --mem 48G
#SBATCH -c 8
#SBATCH --array=10
#SBATCH --time 18:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa
module load samtools

for FILE in `cat x${SLURM_ARRAY_TASK_ID}`; do
	bwa mem -t 8 ./refgenome/ref_genome.fa ./batch1/batch1_demultiplexed/${FILE}.1.fq ./batch1/batch1_demultiplexed/${FILE}.2.fq | samtools view -F 0x04 -f 0x02 -q 20 -h -b -@8 -o ./stacks_align/${FILE}.bam
done
