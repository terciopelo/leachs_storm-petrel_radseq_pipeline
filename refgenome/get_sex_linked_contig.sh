#!/bin/bash
#SBATCH --job-name=align_reference_to_kittiwake
#SBATCH --mem 40G
#SBATCH -c 1
#SBATCH --array 0
#SBATCH --time 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load  minimap2/2.28

#NOTE: minimap2 command must be run in directory that either contains the reference genomes themselves 
##or else symlinks to the reference genomes
#NOTE: individual query contigs are split up by minimap2 for alignment; current filtering strategy throws out contig if one segment 
##aligns to a sex chromosome
#NOTE: Using a conservative mapq score of 13 (corresponding to a 5% chance of incorrect alignment) in order to say 
##particular contig segments do indeed align to a sex chromosome

chromosomal_refGenomeList=("kittiwake_ref.fa")
chromosomal_speciesNames=("rissa_tridactyla")
z_chromosome_names=("NC_071497.1")
w_chromosome_names=("NC_071496.1")

scaffold_refGenomeList=("ref_genome.fa")
scaffold_speciesNames=("pagophila_eburnea")

minimap2 -ax asm5 "${chromosomal_refGenomeList[$SLURM_ARRAY_TASK_ID]}"  "${scaffold_refGenomeList[$SLURM_ARRAY_TASK_ID]}" > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam

awk -v Z="${z_chromosome_names[$SLURM_ARRAY_TASK_ID]}" '($3 == Z && $5 >= 20) {print $1}' "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam | sort -u > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_z_contigs.txt

awk -v W="${w_chromosome_names[$SLURM_ARRAY_TASK_ID]}" '($3 == W && $5 >= 20) {print $1}' "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam | sort -u > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_w_contigs.txt
