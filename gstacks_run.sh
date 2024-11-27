#!/bin/bash
#SBATCH --job-name=gstacks_lounder
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --mem 128G
#SBATCH -c 8
#SBATCH --time 72:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load stacks

gstacks -I ./bams_both_batches -O ./stacks_both --min-mapq 20 -M popmap_w_batches.txt  -t 8 
