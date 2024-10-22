#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH -t 1-0:0:0
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --job-name=ParseFastQ_DBR01
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020 python/2.7.18

python ./ParseDBR_ddRAD/ParseFastQ.py -r Heather_s_Final_Library_2_S1_L003_R1_001.fastq.gz \
         -R Heather_s_Final_Library_2_S1_L003_R2_001.fastq.gz \
         -i ATCACG -e AATT -n ./batch2/batch2_DBR01_R1_output.fastq \
         -N ./batch2/batch2_DBR01_R2_output.fastq \
         --drop ./batch2/batch2_dropped_DBR01.txt -Z -l 0