#!/bin/bash
#SBATCH --job-name=bwa_align_all
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --array=10-17
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH --time 10:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa
module load bcftools
module load samtools
# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST=x${SLURM_ARRAY_TASK_ID}

echo $SAMPLELIST

# loop through sample list
for SAMPLE in `cat $SAMPLELIST`; do
	bwa mem -t 4 refgenome/ref_genome.fa ${SAMPLE}_trimmed.1.fq.gz ${SAMPLE}_trimmed.2.fq.gz | samtools view -b | samtools sort --threads 4 > ${SAMPLE}.bam
done
