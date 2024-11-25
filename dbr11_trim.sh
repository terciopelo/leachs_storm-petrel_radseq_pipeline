#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 3-0:0:0
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=ParseFastQ_DBR11_batch1
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020 python/2.7.18

python ~/ParseDBR_ddRAD/ParseFastQ.py -r Heather_s_Final_Library_S1_L001_R1_001.fastq.gz \
         -R Heather_s_Final_Library_S1_L001_R2_001.fastq.gz \
         -i GGCTAC -e AATT -n ./batch1/batch1_DBR11_R1_output.fastq \
         -N ./batch1/batch1_DBR11_R2_output.fastq \
         --drop ./batch1/batch1_dropped_DBR11.txt -Z -l 3
