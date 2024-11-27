#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=demultiplex_radtags_batch1
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=100G 
#SBATCH -t 48:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020

module load stacks/2.64

process_radtags -1 batch1_R1.fastq \
        -2 batch1_R2.fastq \
        -o ./batch1_demultiplexed \
        -b batch1_barcodes.txt -P -c -r -q --inline-inline \
        --renz_1 sbfI --renz_2 mluCI
