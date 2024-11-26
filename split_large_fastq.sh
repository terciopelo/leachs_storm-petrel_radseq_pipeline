#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 8:00:00
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=split_large_fastq
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# need to first unzip files for processing
# sub in your file names as applicable

gunzip Heather_s_Final_Library_2_S1_L003_R1_001.fastq.gz
gunzip Heather_s_Final_Library_2_S1_L003_R2_001.fastq.gz

# files are the same length, so only need to get the length of one
AVAR=`wc -l < Heather_s_Final_Library_2_S1_L003_R1_001.fastq`

# number of files to generate; can switch if need be
BVAR=8

# get number of lines to split by
LINE_VAR=`echo $((AVAR / BVAR))`

# then split both read files. For later convenience (i.e., when using Slurm task IDs for parallelization), start at 10 and use numeric suffixes
split -l $LINE_VAR --numeric-suffixes=10  --additional-suffix=_R1.fastq Heather_s_Final_Library_2_S1_L003_R1_001.fastq

split -l $LINE_VAR --numeric-suffixes=10  --additional-suffix=_R2.fastq Heather_s_Final_Library_2_S1_L003_R2_001.fastq
