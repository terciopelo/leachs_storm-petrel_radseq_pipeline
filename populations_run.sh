#!/bin/bash
#SBATCH --job-name=populations_lounder
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=11ckb5@queensu.ca
#SBATCH --mem 64G
#SBATCH -c 8
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load stacks

populations -P ./stacks_both -O ./stacks_prelim_test -r 0.75 -R 0.5 --min-maf 0.05 -M popmap_w_batches.txt -t 8 --radpainter --vcf --genepop --structure --plink --phylip --treemix

#module load bwa
#module load bcftools
#module load samtools
# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
#SAMPLELIST=x${SLURM_ARRAY_TASK_ID}

#echo $SAMPLELIST

# loop through sample list
#for SAMPLE in `cat $SAMPLELIST`; do

#	bwa mem -t 4 refgenome/ref_genome.fa ${SAMPLE}_trimmed.1.fq.gz ${SAMPLE}_trimmed.2.fq.gz | samtools view -b | samtools sort --threads 4 > ${SAMPLE}.bam
#done


# filters/analyses done...
# --detect_adapter_for_pe = automatically try to detect adapters in addition to overlap analysis
# --dedup = deduplicate sequences
# --adapter_sequence (and adapter_sequence_r2) = check for sequences specified by Genome Quebec
# --cut_tail = run sliding window trimming from tail of sequence to head to remove low quality area
# -y = low complexity filter (default = 30%)
# -p  = over-representation analysis
# -g = remove polyG tails
# -w = threads 
# -i = input r1
# -I = input r2
# -o = output r1
# -O = output r2
# -h = path for html report
