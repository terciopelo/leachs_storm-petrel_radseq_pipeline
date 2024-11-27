#!/bin/bash
#SBATCH --job-name=populations_lounder
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --mem 64G
#SBATCH -c 8
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load stacks

populations -P ./stacks_both -O ./stacks_prelim_test -r 0.75 -R 0.5 --min-maf 0.05 -M popmap_w_batches.txt -t 8 --radpainter --vcf --genepop --structure --plink --phylip --treemix
