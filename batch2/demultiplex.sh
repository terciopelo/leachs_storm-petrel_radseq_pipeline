#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --job-name=demultiplex_radtag_batch2
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=100G 
#SBATCH -t 16:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020

module load stacks/2.64

process_radtags -1 batch2_all_R1.fastq \
        -2 batch2_all_R2.fastq \
        -o ./batch2_demultiplexed \
        -b batch2_barcodes.txt -P -c -r -q --inline-inline \
        --renz_1 sbfI --renz_2 mluCI
