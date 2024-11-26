#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 8:00:00
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=combine_fastqs_batch2
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e


########Combine all DBRs###########

# combine DBR01s
cat batch2_DBR01_x10_R1_output.fastq batch2_DBR01_x11_R1_output.fastq batch2_DBR01_x12_R1_output.fastq batch2_DBR01_x13_R1_output.fastq batch2_DBR01_x14_R1_output.fastq batch2_DBR01_x15_R1_output.fastq batch2_DBR01_x16_R1_output.fastq batch2_DBR01_x17_R1_output.fastq > batch2_DBR01_all_R1.fastq

cat batch2_DBR01_x10_R2_output.fastq batch2_DBR01_x11_R2_output.fastq batch2_DBR01_x12_R2_output.fastq batch2_DBR01_x13_R2_output.fastq batch2_DBR01_x14_R2_output.fastq batch2_DBR01_x15_R2_output.fastq batch2_DBR01_x16_R2_output.fastq batch2_DBR01_x17_R2_output.fastq > batch2_DBR01_all_R2.fastq

# combine DBR08s
cat batch2_DBR08_x10_R1_output.fastq batch2_DBR08_x11_R1_output.fastq batch2_DBR08_x12_R1_output.fastq batch2_DBR08_x13_R1_output.fastq batch2_DBR08_x14_R1_output.fastq batch2_DBR08_x15_R1_output.fastq batch2_DBR08_x16_R1_output.fastq batch2_DBR08_x17_R1_output.fastq > batch2_DBR08_all_R1.fastq

cat batch2_DBR08_x10_R2_output.fastq batch2_DBR08_x11_R2_output.fastq batch2_DBR08_x12_R2_output.fastq batch2_DBR08_x13_R2_output.fastq batch2_DBR08_x14_R2_output.fastq batch2_DBR08_x15_R2_output.fastq batch2_DBR08_x16_R2_output.fastq batch2_DBR08_x17_R2_output.fastq > batch2_DBR08_all_R2.fastq

# combine DBR10s
cat batch2_DBR10_x10_R1_output.fastq batch2_DBR10_x11_R1_output.fastq batch2_DBR10_x12_R1_output.fastq batch2_DBR10_x13_R1_output.fastq batch2_DBR10_x14_R1_output.fastq batch2_DBR10_x15_R1_output.fastq batch2_DBR10_x16_R1_output.fastq batch2_DBR10_x17_R1_output.fastq > batch2_DBR10_all_R1.fastq

cat batch2_DBR10_x10_R2_output.fastq batch2_DBR10_x11_R2_output.fastq batch2_DBR10_x12_R2_output.fastq batch2_DBR10_x13_R2_output.fastq batch2_DBR10_x14_R2_output.fastq batch2_DBR10_x15_R2_output.fastq batch2_DBR10_x16_R2_output.fastq batch2_DBR10_x17_R2_output.fastq > batch2_DBR10_all_R2.fastq
###################################

#######Combine into omnibus files##

cat batch2_DBR01_all_R1.fastq batch2_DBR08_all_R1.fastq batch2_DBR10_all_R1.fastq > batch2_all_R1.fastq

cat batch2_DBR01_all_R2.fastq batch2_DBR08_all_R2.fastq batch2_DBR10_all_R2.fastq > batch2_all_R2.fast
