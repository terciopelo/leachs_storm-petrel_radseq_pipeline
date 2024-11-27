#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=index_ref_genome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --mem 48G
#SBATCH -c 1
#SBATCH --time 2:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa

bwa index ref_genome.fa
