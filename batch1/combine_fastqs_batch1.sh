#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 8:00:00
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=combine_fastqs_batch1
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# combine all R1s
cat batch1_DBR01_R1_output.fastq batch1_DBR08_R1_output.fastq batch1_DBR10_R1_output.fastq batch1_DBR11_R1_output.fastq > batch1_R1.fastq

# combine all R2s
cat batch1_DBR01_R2_output.fastq batch1_DBR08_R2_output.fastq batch1_DBR10_R2_output.fastq batch1_DBR11_R2_output.fastq > batch1_R2.fastq
