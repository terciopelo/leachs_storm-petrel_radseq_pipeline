## Leach's storm-petrel ddRadSeq pipeline
Prepared by Chris Boccia, Nov. 2024
Most recent update: Jan. 24, 2026


General notes: 
* Though this walkthrough is specific to the Lounder Leach's storm-petrel data set, I've included some information on how to generalize it for other RadSeq data sets.
* For all cluster scripts--I've asked the SLURM engine to send notification emails about job start times, end times, and errors. Since I don't personally want to receive emails about others' jobs, I've replaced all copies of my email address with "your_email". If you would like these notifications, add your email address in place of this.
* This code details how to generate demultiplexed sequences from raw fastq output. If your sequences have already been demultiplexed, skip to Step 4, which details adapter trimming and quality control.

<br>

---

#### Data and code availability

* If need be, the original raw multiplexed fastq files can be made available
  * For Friesen lab members: raw data is available on the lab backup drives
* gstacks data (cleaned, aligned, and called) will be made available following the publication of the data set.
* all SLURM scripts and relevant primer/popmap files can be cloned from this Github repository: https://github.com/terciopelo/leachs_storm-petrel_radseq_pipeline.git


#### Some notes on the data set...

* LESP extraction data set was assembled by Friesen Lab MSc student Heather Lounder
* Library prep for RadSeq was done by Zhengxin Sun in-house at Queen's University
* Degenerate primers were used, following a previous protocol. I've included the explanation files and primer barcode/sample matches that I received from Zhengxin (in the primers folder of the GitHub repo)
* Since this data set was intended originally to be used to test out a new genetic assignment program, some samples overlap between Heather's two batches of data--these can also be used to check/correct for batch effect differences between the two RadSeq runs
* The primer/sample correspondence tables I received had a number of small typos--these have been corrected. However, from Zhengxin's library prep Excel sheets, it appears one sample was submitted in duplicate *within* the first RadSeq batch. Since I realized this after running the demultiplexing script, I believe the current file for the duplicate individual will correspond to sequences linked to its second set of tags (since the first will have been overwritten). I don't think there's any particular reason for having two samples from the same individual from the same RadSeq run, so I haven't bothered going back to re-acquire the 'missing' duplicate sequence.

#### The data set consists of:

* 328 sequences from 308 individuals
  * Samples were sequenced in two batches
  * 20 individuals overlap between batches 1 and 2:


| Sample | Colony
| ---- | ----|
| 1401-47754 | Country |
| 1401-81363 | Kent |
| 1701-14998 | Bon Portage |
| 1701-15029 | Bon Portage |
| 1701-19735 | Country |
| 2461-02977 | Gull |
| 2461-02981 | Gull |
| 2461-07335 | Baccalieu |
| 2461-07387 | Baccalieu |
| 2461-12635 | Baccalieu |
| 2461-12636 | Baccalieu |
| 2461-12657 | Middle Lawn |
| 2461-12669 | Middle Lawn |
| 39334 | Corossol |
| 39338 | Corossol |
| 7901-11063 | Kent|
| 7901-12044 | Gull |
| 8517724 | Hernyken |
| P27 | Green |
| P30 | Green |

Colony breakdown

| Colony | Number of samples
| ----   | ----- |
| Baccalieu | 32 |
| Baja | 3 |
| Bon Portage | 19 |
| Corossol | 14 |
| Country | 21 |
| Green | 19 |
| Gull | 42 |
| Hernyken | 23 |
| Iceland | 15 |
| Kent | 20 |
| Middle Lawn | 16 |
| Unknown | 84 |

#### Basic cluster setup

Log on to your preferred Digital Research Alliance of Canada cluster
```
ssh <username>@beluga.computecanada.ca
```

<br>

Navigate to your scratch folder (large amount of temporary space you can use for large projects such as this)
```
cd scratch
```

<br>

Clone this GitHub repo

```
git clone https://github.com/terciopelo/leachs_storm-petrel_radseq_pipeline.git
```
<br>

Go to wherever your files are backed up/available and run this to transfer them (or use the Globus transfer service)

Sub in your username, and the names of your sequence files if using this code for a different system

```
# batch 1
scp Heather_s_Final_Library_2_S1_L003_R1_001.fastq.gz <your_username>@beluga.computecanada.ca:~/scratch/leachs_storm-petrel_radseq_pipeline.git/.
scp Heather_s_Final_Library_S1_L001_R2_001.fastq.gz <your_username>@beluga.computecanada.ca:~/scratch/leachs_storm-petrel_radseq_pipeline.git/.

# batch 2
scp Heather_s_Final_Library_2_S1_L003_R1_001.fastq.gz <your_username>@beluga.computecanada.ca:~/scratch/leachs_storm-petrel_radseq_pipeline.git/.
scp Heather_s_Final_Library_2_S1_L003_R2_001.fastq.gz <your_username>@beluga.computecanada.ca:~/scratch/leachs_storm-petrel_radseq_pipeline.git/.
```

### Step 1: Deduplication using degenerate primers

* Use custom Python code developed by Evelyn Jensen
  * https://github.com/Eljensen/ParseDBR_ddRAD
  * For more information, see the files in the primer_info folder within this repo

* On Digital Alliance of Canada cluster (or similar)...


```
# Navigate to your home directory to simplify later calls to this package
cd ~
# clone repo
git clone https://github.com/Eljensen/ParseDBR_ddRAD
```

* Four DBRs were used for Batch 1, three for batch 2 (see table for script parameters)

| DBR | -l | -i|
|----| ----| ----|
|DBR01 | 0 | ATCACG |
|DBR08 | 1 | ACTTGA|
|DBR10 | 2 | TAGCTT |
|DBR11 | 3 | GGCTAC |

----

| Batch | DBRs used
| ----| ----|
| Batch 1 | DBR01, DBR08, DBR10, DBR11 |
| Batch 2 | DBR01, DBR08, DBR10 |

* Lounder batch 1 was small enough that it can be processed without additional parallelization
* Lounder batch 2 was larger, and, on a single core, will not run within the 7 day job duration maximum on most Compute Canada clusters
* Instead, divide up the fastq reads into chunks (I found dividing the reads into 8 chunks worked well for batch 2)

#### Deduplicate batch 1

I generated 4 separate DBR script copies (one for each DBR) that will run for 3-4 days each

Only one copy is included below, since all of them are functionally identical aside from DBR parameters and output file names

<br>

This is file dbr01_trim.sh
```
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 4-0:0:0
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=ParseFastQ_DBR01_batch1
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020 python/2.7.18

python ~/ParseDBR_ddRAD/ParseFastQ.py -r Heather_s_Final_Library_S1_L001_R1_001.fastq.gz \
         -R Heather_s_Final_Library_S1_L001_R2_001.fastq.gz \
         -i ATCACG -e AATT -n /batch1/batch1_DBR01_R1_output.fastq \
         -N ./batch1/batch1_DBR01_R2_output.fastq \
         --drop ./batch1/batch1_dropped_DBR01.txt -Z -l 0

```
#### Notes on above file

* If using on a different dataset, you need to sub in your fastq file names (R1 and R2) into each unparallelized DBR script (dbr01_trim.sh, dbr08_trim.sh, dbr10_trim.sh, drb11_trim.sh)

* The Jensen ParseFastQ.py script requires an older version of python (and also an older version of the standard cluster environment--see the module load statement in the script)

* The -Z flag tells the program not to zip the files after processing. This uses up more space, but speeds up the program, and also speeds up later processing to a degree, since we need to combine all the DBR files (though, if you do zip them, this can be done with zcat)

<br>

#### Run DBR deduplication

Enter the cloned GitHub repository
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
```

<br>

Launch each DBR job!

```
sbatch dbr01_trim.sh
sbatch dbr08_trim.sh
sbatch dbr10_trim.sh
sbatch dbr11_trim.sh
```

These will output all deduplicated sequences to your batch 1 folder. Each script will generate an R1 and R2 file containing all sequences for that DBR, along with outputting all dropped reads to a "_dropped" file in case those are of interest.

Since batch 2 had a larger number of sequences, we need to...

#### Split large batch 2 fastq files into smaller chunks for parallel processing

* Note: fastq files are arranged in sets of four lines, so you can easily split them apart for parallel processing. The correspondence between lines in R1 and R2 should be exact

* When dividing, you need to ensure your divisor is itself divisible by 4, and even divisors are preferable (e.g., 8, 16)

* If this operation times out in a salloc, I've made an automated SLURM script (split_large_fastq.sh)

```
# (probably run the code below in an interactive allocation or slurm script to avoid overusing the login nodes)
salloc --mem 40G --time 2:00:00

# need to first unzip files for processing
gunzip Heather_s_Final_Library_2_S1_L003_R1_001.fastq.gz
gunzip Heather_s_Final_Library_2_S1_L003_R2_001.fastq.gz

# determine length of files 
wc -l Heather_s_Final_Library_2_S1_L003_R1_001.fastq

# divide the number obtained above by 8 (generate 8 files per read direction for parallelization)

# then split both read files. For later convenience (i.e., when using Slurm task IDs for parallelization), start at 10 and use numeric suffixes
split -l 1631171478 --numeric-suffixes=10  --additional-suffix=_R1.fastq Heather_s_Final_Library_2_S1_L003_R1_001.fastq

split -l 1631171478 --numeric-suffixes=10  --additional-suffix=_R2.fastq Heather_s_Final_Library_2_S1_L003_R2_001.fastq
```
Notes:
* If using this on a different system, sub the number you got when diving by 8 in place of "1631171478"
* If the result of wc -l is less than 1631171478, you likely don't need to split the fastq files, and can refer to the code for batch 1
* This code will generate 8 files labelled x10-x17_R{1/2}.fastq

<br>

Deduplicate batch 2

* Example included below. 

```
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 4-0:0:0
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=ParseFastQ_DBR01_batch2array
#SBATCH --array=10-17
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020 python/2.7.18

python ~/ParseDBR_ddRAD/ParseFastQ.py -r x${SLURM_ARRAY_TASK_ID}_R1.fastq  \
         -R x${SLURM_ARRAY_TASK_ID}_R2.fastq \
         -i ATCACG -e AATT -n ./batch2/batch2_DBR01_x${SLURM_ARRAY_TASK_ID}_R1_output.fastq \
         -N ./batch2/batch2_DBR01_x${SLURM_ARRAY_TASK_ID}_R2_output.fastq \
         --drop ./batch2/batch2_dropped_DBR01_${SLURM_ARRAY_TASK_ID}.txt -Z -l 0

```
Notes
* This script will launch a separate job for each of the 8 different paired chunks of the batch 2 fastq file (i.e., x##_R1.fastq and x##_R2.fastq)
* The --array=10-17 slurm flag causes the SLURM_ARRAY_TASK_ID variable to take on the numbers between 10 and 17, which match our split files (x10-17)
* The script similarly outputs two files tagged with both the DBR and chunk number

<br>

Run DBR array scripts (only three this time)

```
sbatch dbr01_trim_array.sh
sbatch dbr08_trim_array.sh
sbatch dbr10_trim_array.sh
```

### Step 2: Combining files

In order to run the demultiplexing procedure, need to have two omnibus read files per batch

<br>

#### Batch 1

Combine all deduplication fastqs into omnibus fastq files
* If this times out, can run using a SLURM script, see (batch1/combine_fastqs_batch1.sh) 
* You can use wildcard flags for this, but I'm always a little concerned the sorting for those might be weird...*Should* be fine though

```
cd batch1

salloc --mem 40G --time 2:00:00

# combine all R1s
cat batch1_DBR01_R1_output.fastq batch1_DBR08_R1_output.fastq batch1_DBR10_R1_output.fastq batch1_DBR11_R1_output.fastq > batch1_R1.fastq

# wildcard version
#cat batch1_DBR*_R1_output.fastq > batch1_R1.fastq

# combine all R2s
cat batch1_DBR01_R2_output.fastq batch1_DBR08_R2_output.fastq batch1_DBR10_R2_output.fastq batch1_DBR11_R2_output.fastq > batch1_R2.fastq

# wildcard version
#cat batch1_DBR*_R2_output.fastq > batch1_R2.fastq
```

#### Batch 2

For batch 2, there's an extra step--we have to combine all the parallelized fastq file chunks for each DBR before we can do our final concatenation

Assuming you're still in your batch1 folder...
* Again, wildcards are probably better here, if you're running this on a different data set
* If running on the Lounder data set, can use (batch2/combine_fastqs_batch2.sh)

```
cd ../batch2

# combine DBR01s
cat batch2_DBR01_x10_R1_output.fastq batch2_DBR01_x11_R1_output.fastq batch2_DBR01_x12_R1_output.fastq batch2_DBR01_x13_R1_output.fastq batch2_DBR01_x14_R1_output.fastq batch2_DBR01_x15_R1_output.fastq batch2_DBR01_x16_R1_output.fastq batch2_DBR01_x17_R1_output.fastq > batch2_DBR01_all_R1.fastq
cat batch2_DBR01_x10_R2_output.fastq batch2_DBR01_x11_R2_output.fastq batch2_DBR01_x12_R2_output.fastq batch2_DBR01_x13_R2_output.fastq batch2_DBR01_x14_R2_output.fastq batch2_DBR01_x15_R2_output.fastq batch2_DBR01_x16_R2_output.fastq batch2_DBR01_x17_R2_output.fastq > batch2_DBR01_all_R2.fastq

# combine DBR08s
cat batch2_DBR08_x10_R1_output.fastq batch2_DBR08_x11_R1_output.fastq batch2_DBR08_x12_R1_output.fastq batch2_DBR08_x13_R1_output.fastq batch2_DBR08_x14_R1_output.fastq batch2_DBR08_x15_R1_output.fastq batch2_DBR08_x16_R1_output.fastq batch2_DBR08_x17_R1_output.fastq > batch2_DBR08_all_R1.fastq
cat batch2_DBR08_x10_R2_output.fastq batch2_DBR08_x11_R2_output.fastq batch2_DBR08_x12_R2_output.fastq batch2_DBR08_x13_R2_output.fastq batch2_DBR08_x14_R2_output.fastq batch2_DBR08_x15_R2_output.fastq batch2_DBR08_x16_R2_output.fastq batch2_DBR08_x17_R2_output.fastq > batch2_DBR08_all_R2.fastq

# combine DBR10s
cat batch2_DBR10_x10_R1_output.fastq batch2_DBR10_x11_R1_output.fastq batch2_DBR10_x12_R1_output.fastq batch2_DBR10_x13_R1_output.fastq batch2_DBR10_x14_R1_output.fastq batch2_DBR10_x15_R1_output.fastq batch2_DBR10_x16_R1_output.fastq batch2_DBR10_x17_R1_output.fastq > batch2_DBR10_all_R1.fastq
cat batch2_DBR10_x10_R2_output.fastq batch2_DBR10_x11_R2_output.fastq batch2_DBR10_x12_R2_output.fastq batch2_DBR10_x13_R2_output.fastq batch2_DBR10_x14_R2_output.fastq batch2_DBR10_x15_R2_output.fastq batch2_DBR10_x16_R2_output.fastq batch2_DBR10_x17_R2_output.fastq > batch2_DBR10_all_R2.fastq

#### Combine these files into one omnibus fastq pair
cat batch2_DBR01_all_R1.fastq batch2_DBR08_all_R1.fastq batch2_DBR10_all_R1.fastq > batch2_all_R1.fastq
cat batch2_DBR01_all_R2.fastq batch2_DBR08_all_R2.fastq batch2_DBR10_all_R2.fastq > batch2_all_R2.fastq
```

### Step 3: Demultiplexing using Stacks

This is done separately for batch 1 and batch 2

<br>

* For this, you will need some of the information mentioned in the "Preparation" section:
  * Your digestion enzymes (see https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php for the full list)
    * For the Lounder data set, these are...
      * sbfI 
      * mluCI
  * Your "barcode file", linking sample IDs to primer barcodes
    * This should be in the following format, the {tab}s replaced by hitting the tab key:
      * adapter1{tab}adapter2{tab}sample_ID
        * e.g., 
        ```
        CTCG  ATCACG  1401-47758
        TGCA  ATCACG  1401-47748
        ...
        ATCGTA  GGCTAC  2461-13042
        CATCGT  GGCTAC  2461-13480
        ```

        * Barcode files are headerless (i.e., don't put the column names at the top)
        * Zhengxin usually  provides this information in a pair of files--one which lists the possibilities for adapter 1, which is dependent on the enzymes used, and a second which lists the adapter used from the aforementioned list in the first column, and the degenerate DBR primer used in the second column. This pairing ensures uniquely barcoded reads
        * If running on a different data set: once you have generated this file, copy it into the 'batch1' folder
        * For convenience, you can name it 'batch1_barcodes.txt', the name I've used for the Lounder batch 1 barcode file. You can also view this file to check if you've correctly formatted your own barcode file

#### Batch 1

To demultiplex all individuals in batch 1, first make an output directory within batch1

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline/batch1
mkdir batch1_demultiplexed
```

Then run the following script:
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=demultiplex_radtags
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

```

Notes
* Once again, like with python2, Stacks/2.64 requires loading StdEnv/2020, an older environment module
* the --inline-inline flag means that barcodes are present on both forward and reverse reads

```
sbatch demultiplex.sh
```

#### Batch 2
Essentially the same as batch 1.

To demultiplex all individuals in batch 2, first make an output directory within batch2

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline/batch2
mkdir batch2_demultiplexed
```

Then run the following script:
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=demultiplex_radtags_batch2
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 48:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020

module load stacks/2.64

process_radtags -1 batch2_all_R1.fastq \
        -2 batch2_all_R2.fastq \
        -o ./batch2_demultiplexed \
        -b batch2_barcodes.txt -P -c -r -q --inline-inline \
        --renz_1 sbfI --renz_2 mluCI

```

```
sbatch demultiplex.sh
```

The output from process_radtags is comprised of pairs of gzipped fastq files, one pair for each sample in the batch. File names are formatted as follows:
* samplename.1.fq.gz 
* samplename.2.fq.gz

### Step 4: Trimming and quality control

Next, it's probably a good idea to run fastp to remove any latent adapter sequences or low quality bases. You may need to do this a couple of times to get the parameters right. I noticed that many of these sequences had repetitive bases at the start and end, so I added some flags to remove these. Fastp removes adapter content, low quality bases, and can accommodate a variety of other quality control scans/adjustments. It also outputs its own version of a 'fastqc' sequence quality report. You can view these for each sequence (look for the .html files that match each sample name).

First, we need to generate some files listing all of our samples and their paths so that we can run fastp in parallel

<br>

Get sample list for batch 1
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
ls -1 batch1/batch1_demultiplexed/*1.fq.gz | sed 's/........$//' > all_batch1_demultiplexed.txt
```
Note: sed removes the last 8 characters, so we get just the sample name, not the file ending, which is helpful.

Get sample list for batch 2
```
ls -1 batch2/batch2_demultiplexed/*1.fq.gz | sed 's/........$//' > all_batch2_demultiplexed.txt
```

Combine the files

```
cat all_batch1_demultiplexed.txt all_batch2_demultiplexed.txt > all_demultiplexed_sample_paths.txt
```

Generate 8 roughly even splits
```
AVAR=`wc -l < all_demultiplexed_sample_paths.txt`

# number of files to generate; can switch if need be
BVAR=8

# get number of lines to split by
LINE_VAR=`echo $((AVAR / BVAR))`

# then split both read files. For later convenience (i.e., when using Slurm task IDs for parallelization), start at 10 and use numeric suffixes
split -l $LINE_VAR --numeric-suffixes=10  all_demultiplexed_sample_paths.txt
```
Note: since LINE_VAR must be an integer, bash uses the div operation, rather than float/double division which permits decimal points. This means that if your "number of files" is not a factor of the number of samples that you have, you will end up with one additional file that contains the remainder of the samples (e.g., if you have 42 samples and want 8 splits, this script will generate 8 files with 5 samples each and a ninth file with 2 samples).

#### Run fastp

Run the following array script. It will generate 8 different copies of the script, one for each of the sample splits done above. Trimmed files will be saved in the same directory as the raw demultiplxed files (batch1/batch1_demultiplexed and batch2/batch2_demultiplexed respectively.

```
#!/bin/bash
#SBATCH --job-name=fastp_lesp_all
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --array=10-17
#SBATCH --mem 20G
#SBATCH -c 4
#SBATCH --time 10:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastp

# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST=x${SLURM_ARRAY_TASK_ID}

echo $SAMPLELIST

# loop through sample list
for SAMPLE in `cat $SAMPLELIST`; do
        fastp --detect_adapter_for_pe -w 4 --dedup -f 8 -F 4 -i ${SAMPLE}.1.fq.gz -I ${SAMPLE}.2.fq.gz -o ${SAMPLE}_trimmed.1.fq.gz -O ${SAMPLE}_trimmed.2.fq.gz -h ${SAMPLE}.html
done

# filters/analyses done...
# --detect_adapter_for_pe = automatically try to detect adapters in addition to overlap analysis
# --dedup = deduplicate sequences
# -f 8 = cut first 8 bases
# -F 4 = cut trailing 4 bases
```

```
sbatch fastp_trim.sh
```

Notes:
* Some of these filters/cleaning flags are likely not necessary if you're running a different data set through this pipeline.
* I would personally try running this without the two -f and -F filters first. Have a look at the resulting html sequence quality reports, and then see if any additional trimming parameters are needed; then re-run the script if necessary

### Step 5: Alignment and sorting

Next, use bwa mem to align sequences to the Leach's storm-petrel reference genome.

First, enter the reference genome folder

```
cd refgenome
```

Next, download the reference genome from NCBI

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/449/065/GCA_030449065.1_OLeu_1.0/GCA_030449065.1_OLeu_1.0_genomic.fna.gz
```

Unzip the file and rename it for convenience
```
gunzip GCA_030449065.1_OLeu_1.0_genomic.fna.gz
mv GCA_030449065.1_OLeu_1.0_genomic.fna ref_genome.fa
```

Index the reference genome with bwa using the following script

```
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
```
```
sbatch index_refgenome.sh
```
Notes:
* If using this code on a different species, find your species' reference genome on NCBI genome (https://www.ncbi.nlm.nih.gov/home/genomes/) and download/find the download link for your reference genome.

#### Run bwa mem

Run the following array script that uses the same parallelization as the trimming script

```
#!/bin/bash
#SBATCH --job-name=bwa_align_all
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --array=10-17
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH --time 10:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa
module load bcftools
module load samtools
# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST=x${SLURM_ARRAY_TASK_ID}

echo $SAMPLELIST

# loop through sample list
for SAMPLE in `cat $SAMPLELIST`; do
        bwa mem -t 4 refgenome/ref_genome.fa ${SAMPLE}_trimmed.1.fq.gz ${SAMPLE}_trimmed.2.fq.gz | samtools view -b | samtools sort --threads 4 > ${SAMPLE}.bam
done
```

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
sbatch align_bwa.sh
```

Notes:
* As with the trimming script, this script will output alignment files to the same directory as the source trimmed fastq files
* This step combines forward and reverse reads, so you get only a single output file in bam format
* The script also sorts each file, which I believe is necessary for Stacks (or at least makes running Stacks faster)

First, since some samples overlap between batches, we need to add a prefix to each bam file

Batch 1
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
cd batch1/batch1_demultiplexed

# for each bam file, add b1_ to the start of the file name
for FILE in `ls -1 *.bam`; do
 mv $FILE b1_${FILE}
done
```

Batch 2
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
cd batch2/batch2_demultiplexed

# for each bam file, add b2_ to the start of the file name
for FILE in `ls -1 *.bam`; do
 mv $FILE b2_${FILE}
done
```

Note:
* This necessitates that your popmap file, which links samples to populations, also has these prefixes. I've added these to the popmap_w_batches.txt file
* Ideally you won't need this step if using this pipeline on other species

Finally, generate a folder to house all bam files, and move them there. Also make an output directory for Stacks

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
mkdir bams_both_batches
mv batch2/batch2_demultiplexed/*.bam bams_both_batches/.
mv batch1/batch1_demultiplexed/*.bam bams_both_batches/.
mkdir stacks_both
```

#### Run gstacks

Run the following script to generate a Stacks output folder, which includes variant calls (identification of variant sites within Rad loci). You can then use this output folder to generate filtered datasets; filtering is very quick and lightweight, so you can sometimes do it without using a SLURM script.

```
#!/bin/bash
#SBATCH --job-name=gstacks_lounder
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --mem 128G
#SBATCH -c 8
#SBATCH --time 72:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load stacks

gstacks -I ./bams_both_batches -O ./stacks_both --min-mapq 20 -M popmap_w_batches.txt  -t 8
```

### Step 7: Run populations to filter your data

Filtering RadSeq data can be tricky. The parameters I use below are preliminary for this data set. For a more extensive treatment of the topic see this Speciation Genomics walkthrough (https://speciationgenomics.github.io/filtering_vcfs/) that provides code for assessing summary stats and filtering accordingly.

#### Depth filtering
* Extremely high depth sites may be mitochondrial, or may be duplicated regions, or potentially even sequening errors, so it's probably a good idea to identify them and remove them. Similarly, low depth sites may not be called with confidence.

* Probably the most direct way of doing this is to run populations preliminarily (since populations unfortunately doesn't have a max depth thresholding option)

Run populations
e.g.,
```
salloc --mem 40G --time 2:00:00 -c 8
module load stacks
populations -P ./stacks_both -O ./stacks_prelim_test -M popmap_w_batches.txt -t 8 --radpainter --vcf --vcf-all --genepop --structure --plink --phylip --treemix
```

After you've run populations, use vcftools to obtain a depth-per-site report, averaged by individual

```
cd stacks_prelim_test
module load vcftools
vcftools --site-mean-depth --vcf populations.all.vcf --out depth_stats
```

Next, go to your gstacks folder, and obtain depth statistics from gstacks logfile
```
grep coverage ../stacks_both/gstacks.log
```

This will yield, for the Leach's storm-petrel data set, for example...
```
effective per-sample coverage: mean=28.9x, stdev=19.7x, min=1.2x, max=90.5x
```

It probably makes sense to set a hard lower cutoff at 8x, and then have the upper threshold be the mean + 2 SD (68.3x)

If interested, you can also open the depth file in R and generate a depth histogram

```
module load r/4.4.0
R

library(tidyverse)
dat = read_delim("depth_stats.ldepth.mean")

num_breaks = max(dat$MEAN_DEPTH,na.rm=T)

# output a PDF plot
pdf("depth_stats_persite_histogram.pdf")
hist(dat$MEAN_DEPTH, breaks = num_breaks)
dev.off()
```

Use this histogram to assess your thresholding choices.

Once you're satisfied with your depth thresholds, run the following command to obtain scaffolds with unusually high or low depth:

```
awk '$3 >= 68.3 || $3 < 8' depth_stats.ldepth.mean > depth_droplist.txt

# since the 'all sites' VCF has a lot of NAs, drop these for downstream convenience
awk '$3 != "-nan"' depth_droplist.txt > depth_droplist_nonan.txt
```

Now, switch to working on the the Stacks catalog to figure out which loci correspond to these scaffold positions...

```
# navigate to your gstacks folder
cd ../stacks_both

# get headers, convert to a bed-file-esque format
# get header lines only
zcat catalog.fa.gz  | awk '(NR % 2) {print}'> catalog_headers.txt

# trim to only relevant fields (chromosome, start position)
cut -d: -f 1,2 catalog_headers.txt > catalog_headers_trimmed.txt

# drop other superfluous characters
sed 's/>//g' -i catalog_headers_trimmed.txt
sed 's/pos=//g' -i catalog_headers_trimmed.txt
sed $'s/:/\t/g' -i catalog_headers_trimmed.txt
sed $'s/ /\t/g' -i catalog_headers_trimmed.txt
```

Get lengths of sequences from the non-header fasta lines
```
# get lengths of loci
zcat catalog.fa.gz  | awk '!(NR % 2) {print length($0)}'> catalog_lengths.txt
```

Get positive or negative strand...
```
cut -d: -f 3 catalog_headers.txt | cut -d " " -f 1 > catalog_headers_strand.txt
```

Combine all of this into one file
```
paste catalog_headers_trimmed.txt catalog_lengths.txt catalog_headers_strand.txt > catalog_loci_lookup.txt
```

Run a quick r-script to convert positive/negative strand start/ends to unidirectional (since we don't need the actual bases/directionality)

```
R
library(tidyverse)
catal = read_delim("./catalog_loci_lookup.txt", col_names=F)

# split catalog into different strand directions
minus  = catal %>% filter(X5 == "-")
plus = catal %>% filter(X5 == "+")

# create start and end variables as appropriate
plus$start = plus$X3
plus$end = plus$X3 + plus$X4
minus$start = minus$X3 - minus$X4
minus$end = minus$X3

# recombine and sort by locus ID
combined = rbind(plus,minus)
combined = combined %>% arrange(X1)

write_delim(combined, delim="\t", col_names=F, file = "catalog_strand_lookup.txt")
```

Finally, generate the depth blacklist!

```
library(tidyverse)

to_rem = read_delim("../stacks_prelim_test/depth_droplist_nonan.txt")
catal = read_delim("./catalog_strand_lookup.txt", col_names=F)
remove_loci= c()
for(i in 1:length(to_rem$CHROM)){
   temp = catal %>% filter(X2 == to_rem$CHROM[i]) %>% filter(X6 <= to_rem$POS[i] & X7 >=to_rem$POS[i])
   remove_loci = c(remove_loci, temp$X1)
}

# since there will be some double-hits...
remove_loci = unique(remove_loci)

r_l = data.frame(rl=remove_loci)
write_delim(r_l, delim="\t", col_names=F, file="depth_blacklist.txt")
```

You can combine this blacklist with the one generated in the section below on removing sex chromosome linked sites via

```
cat depth_blacklist.txt blacklist.txt > ../refgenome/overall_blacklist.txt
```

I've arbitrarily decided to put the blacklist in the refgenome folder.

#### Identifying scaffolds on sex chromosomes
* This is still a bit experimental. Thanks to Spencer for the original code!
* This script will align your reference genome against a reference genome of closely-related species that contains sex chromosomes
* If your reference genome already contains sex chromosomes, there's no need to run the aligment script--you can just take the names of the sex chromosomes in your alignment and skip to blacklist generation
* Should update this to include mitochondrial loci as well?

First, obtain a close-relative reference genome from NCBI that is chromosome-level (you can use the "Taxonomy" feature to help you find this)

Enter your reference genome folder and download/unzip the genome
e.g., for the kittiwake chromosome-level reference genome...

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline/refgenome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/500/815/GCF_028500815.1_bRisTri1.patW.cur.20221130/GCF_028500815.1_bRisTri1.patW.cur.20221130_genomic.fna.gz

gunzip GCF_028500815.1_bRisTri1.patW.cur.20221130_genomic.fna.gz
mv GCF_028500815.1_bRisTri1.patW.cur.20221130_genomic.fna kittiwake_ref.fa
```

Next, adjust find the names of the sex chromosomes (can do this on NCBI by clicking directly on them on the genome page). Add them to the "get_sex_linked_contig.sh" script, and adjust the other variables in the script according to your data set.
```
chromosomal_refGenomeList=("kittiwake_ref.fa") #this is the genome to align against (the chromosome-level one)
chromosomal_speciesNames=("rissa_tridactyla") # name of the species for the chromosome-level reference
z_chromosome_names=("NC_071497.1") # z chromosome ID
w_chromosome_names=("NC_071496.1") # w chromosome ID

scaffold_refGenomeList=("ref_genome.fa") # your focal species' reference genome (generally no need to change)
scaffold_speciesNames=("pagophila_eburnea") # name of your focal species
```

Next, run the script:

```
sbatch get_sex_linked_contig.sh
```

You will use then use the resulting contig files and your gstacks catalog to identify "Stacks loci" to remove. The correspondence between Stacks loci numbers and reference genome scaffolds is present in the catalog.fa.gz file output by gstacks.

e.g., for the kittiwake/ivory gull sex chromosome alignment...
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline/refgenome

cat *_w_contigs.txt *_z_contigs.txt > all_sex_linked_contigs.txt

touch blacklist.txt
for CONTIG in `cat all_sex_linked_contigs.txt`; do
 #change path to catalog file to reflect what you named your gstacks output folder
 zgrep ${CONTIG} ../stacks_original/catalog.fa.gz | cut -f 1 -d " "| sed 's/>//g'  > blacklist.txt
done
sort -o blacklist.txt blacklist.txt
```
This file can now be used to remove sex-linked loci from your populations output.


First, make an output directory

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
mkdir stacks_prelim_test
```

Run script

```
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

populations -P ./stacks_both -O ./stacks_prelim_test -r 0.75 -R 0.5 --blacklist ./refgenome/overall_blacklist.txt --min-maf 0.05 -M popmap_w_batches.txt -t 8 --radpainter --vcf --vcf-all --genepop --structure --plink --phylip --treemix
```

```
sbatch populations_run.sh
```
Notes:
* This script will generate input formats for a wide variety of downstream programs. VCF files are the most widely used, but radpainter, structure, plink, and genepop are also commonly used. Phylip and treemix are mainly used for explicit phylogenetic reconstruction and comparative methods.

* Note: for Stairway plot analyses, omit the --min-maf 0.05 flag and parameter

### Steps before downstream processing

At this point, most folks use easySFS to generate a site frequency spectrum, and plink or similar to do PCA

Once you have determined the appropriate population structure to use (i.e., if you are planning to combine populations that appear to be one metapopulation, etc.), generate a new popmap file that captures these changes. This walkthrough will imagine that you have named such a file "leachs_onepop.popmap". For the Leach's data set, I have currently done this by removing all unknowns, and the samples that overlapped between batch 1 and 2 (removed the ones from batch 1), and setting all population names to LEACH. If dropping populations, it's a good idea to re-run populations with the new popmap since you'll likely be able to retain more sites.

### Step 8: Generating a site frequency spectrum using easySFS

Clone easySFS GitHub repo
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
git clone https://github.com/isaacovercast/easySFS.git
```

Set up a salloc and load some modules
* easySFS usually runs very quickly, but should still use an interactive allocation
* Need to load 'scipy-stack' in order to use numPy, other Python libraries
```
salloc --mem 40G --time 2:00:00
module load python
module load scipy-stack/2024a
```

Now, run easySFS using your vcf from populations. 

First, can consider 'projecting' to a different sample of haploid sequences, Maximizing number of segregating sites is ideal, but might have the unintended consequence of adjusting the appropriate '0' bin ('L' parameter).
```
python easySFS/easySFS.py -i ./stacks_prelim_test/populations.snps.vcf -p leachs_onepop.popmap -a --preview
```

Sample output (from Leach's data):
```
LEACH
(2, 8366)       (3, 12515)      (4, 15020)      (5, 16661)      (6, 17856)      (7, 18734)      (8, 19432)      (9, 19971)      (10, 20427)     (11, 20793)     (12, 21113) (13, 21369)      (14, 21604)     (15, 21791)     (16, 21969)     (17, 22106)     (18, 22244)     (19, 22355)     (20, 22464)     (21, 22551)     (22, 22639)     (23, 22707) (24, 22778)      (25, 22828)     (26, 22886)     (27, 22929)     (28, 22977)     (29, 23003)     (30, 23043)     (31, 23065)     (32, 23098)     (33, 23109)     (34, 23137) (35, 23140)      (36, 23163)     (37, 23167)     (38, 23187)     (39, 23178)     (40, 23195)     (41, 23190)     (42, 23204)     (43, 23192)     (44, 23204)     (45, 23201) (46, 23211)      (47, 23195)     (48, 23204)     (49, 23189)     (50, 23196)     (51, 23178)     (52, 23184)     (53, 23155)     (54, 23161)     (55, 23141)     (56, 23145) (57, 23114)      (58, 23118)     (59, 23086)     (60, 23090)     (61, 23062)     (62, 23065)     (63, 23035)     (64, 23037)     (65, 23007)     (66, 23010)     (67, 22972) (68, 22973)      (69, 22938)     (70, 22940)     (71, 22902)     (72, 22904)     (73, 22875)     (74, 22876)     (75, 22832)     (76, 22833)     (77, 22791)     (78, 22792) (79, 22746)      (80, 22747)     (81, 22694)     (82, 22694)     (83, 22644)     (84, 22644)     (85, 22614)     (86, 22614)     (87, 22565)     (88, 22565)     (89, 22517) (90, 22517)      (91, 22473)     (92, 22474)     (93, 22425)     (94, 22425)     (95, 22375)     (96, 22375)     (97, 22321)     (98, 22321)     (99, 22274)     (100, 22275)(101, 22223)     (102, 22223)    (103, 22175)    (104, 22175)    (105, 22133)    (106, 22133)    (107, 22082)    (108, 22082)    (109, 22029)    (110, 22029)    (111, 21970)(112, 21970)     (113, 21900)    (114, 21900)    (115, 21846)    (116, 21846)    (117, 21786)    (118, 21786)    (119, 21734)    (120, 21734)    (121, 21675)    (122, 21675)(123, 21619)     (124, 21620)    (125, 21562)    (126, 21562)    (127, 21514)    (128, 21514)    (129, 21443)    (130, 21443)    (131, 21376)    (132, 21376)    (133, 21293)(134, 21293)     (135, 21242)    (136, 21242)    (137, 21177)    (138, 21177)    (139, 21092)    (140, 21092)    (141, 21041)    (142, 21041)    (143, 20990)    (144, 20990)(145, 20920)     (146, 20920)    (147, 20873)    (148, 20873)    (149, 20816)    (150, 20816)    (151, 20751)    (152, 20751)    (153, 20684)    (154, 20684)    (155, 20620)(156, 20620)     (157, 20548)    (158, 20548)    (159, 20478)    (160, 20478)    (161, 20417)    (162, 20417)    (163, 20339)    (164, 20339)    (165, 20273)    (166, 20273)(167, 20212)     (168, 20212)    (169, 20139)    (170, 20139)    (171, 20059)    (172, 20059)    (173, 19978)    (174, 19978)    (175, 19901)    (176, 19901)    (177, 19810)(178, 19810)     (179, 19734)    (180, 19734)    (181, 19667)    (182, 19667)    (183, 19597)    (184, 19597)    (185, 19515)    (186, 19515)    (187, 19460)    (188, 19460)(189, 19393)     (190, 19393)    (191, 19313)    (192, 19313)    (193, 19241)    (194, 19241)    (195, 19166)    (196, 19166)    (197, 19100)    (198, 19100)    (199, 19027)(200, 19027)     (201, 18935)    (202, 18935)    (203, 18856)    (204, 18856)    (205, 18793)    (206, 18793)    (207, 18712)    (208, 18712)    (209, 18646)    (210, 18646)(211, 18559)     (212, 18559)    (213, 18482)    (214, 18482)    (215, 18412)    (216, 18412)    (217, 18333)    (218, 18333)    (219, 18248)    (220, 18248)    (221, 18158)(222, 18158)     (223, 18073)    (224, 18073)    (225, 17993)    (226, 17993)    (227, 17909)    (228, 17909)    (229, 17826)    (230, 17826)    (231, 17760)    (232, 17760)(233, 17664)     (234, 17664)    (235, 17600)    (236, 17600)    (237, 17509)    (238, 17509)    (239, 17428)    (240, 17428)    (241, 17334)    (242, 17334)    (243, 17240)(244, 17240)     (245, 17145)    (246, 17145)    (247, 17055)    (248, 17055)    (249, 16961)    (250, 16961)    (251, 16870)    (252, 16870)    (253, 16790)    (254, 16790)(255, 16708)     (256, 16708)    (257, 16623)    (258, 16623)    (259, 16520)    (260, 16520)    (261, 16438)    (262, 16438)    (263, 16349)    (264, 16349)    (265, 16267)(266, 16267)     (267, 16164)    (268, 16164)    (269, 16070)    (270, 16070)    (271, 15981)    (272, 15981)    (273, 15891)    (274, 15891)    (275, 15797)    (276, 15797)(277, 15692)     (278, 15692)    (279, 15590)    (280, 15590)    (281, 15497)    (282, 15497)    (283, 15401)    (284, 15401)    (285, 15308)    (286, 15308)    (287, 15193)(288, 15193)     (289, 15104)    (290, 15104)    (291, 15027)    (292, 15027)    (293, 14921)    (294, 14921)    (295, 14828)    (296, 14828)    (297, 14717)    (298, 14717)(299, 14611)     (300, 14611)    (301, 14512)    (302, 14512)    (303, 14413)    (304, 14413)    (305, 14326)    (306, 14326)    (307, 14231)    (308, 14231)    (309, 14149)(310, 14149)     (311, 14062)    (312, 14062)    (313, 13939)    (314, 13939)    (315, 13825)    (316, 13825)    (317, 13704)    (318, 13704)    (319, 13609)    (320, 13609)(321, 13522)     (322, 13522)    (323, 13444)    (324, 13444)    (325, 13336)    (326, 13336)    (327, 13232)    (328, 13232)    (329, 13127)    (330, 13127)    (331, 13032)(332, 13032)     (333, 12916)    (334, 12916)    (335, 12809)    (336, 12809)    (337, 12652)    (338, 12652)    (339, 12497)    (340, 12497)    (341, 12354)    (342, 12354)(343, 12226)     (344, 12226)    (345, 12097)    (346, 12097)    (347, 11946)    (348, 11946)    (349, 11820)    (350, 11820)    (351, 11687)    (352, 11687)    (353, 11546)(354, 11546)     (355, 11401)    (356, 11401)    (357, 11256)    (358, 11256)    (359, 11106)    (360, 11106)    (361, 10976)    (362, 10976)    (363, 10828)    (364, 10828)(365, 10693)     (366, 10693)    (367, 10549)    (368, 10549)    (369, 10400)    (370, 10400)    (371, 10266)    (372, 10266)    (373, 10117)    (374, 10117)    (375, 9940) (376, 9940)      (377, 9800)     (378, 9800)     (379, 9645)     (380, 9645)     (381, 9481)     (382, 9481)     (383, 9357)     (384, 9357)     (385, 9207)     (386, 9207) (387, 9042)      (388, 9042)     (389, 8887)     (390, 8887)     (391, 8742)     (392, 8742)     (393, 8596)     (394, 8596)     (395, 8457)     (396, 8457)     (397, 8315) (398, 8315)      (399, 8202)     (400, 8202)     (401, 8091)     (402, 8091)     (403, 7998)     (404, 7998)     (405, 7908)     (406, 7908)     (407, 7830)     (408, 7830) (409, 7738)      (410, 7738)     (411, 7603)     (412, 7603)     (413, 7493)     (414, 7493)     (415, 7366)     (416, 7366)     (417, 7235)     (418, 7235)     (419, 7093) (420, 7093)      (421, 6931)     (422, 6931)     (423, 6740)     (424, 6740)     (425, 6532)     (426, 6532)     (427, 6303)     (428, 6303)     (429, 6039)     (430, 6039) (431, 5683)      (432, 5683)     (433, 5248)     (434, 5248)     (435, 4755)     (436, 4755)     (437, 4160)     (438, 4160)     (439, 3490)     (440, 3490)     (441, 2652) (442, 2652)      (443, 1717)     (444, 1717)     (445, 742)      (446, 742)      (447, 96)       (448, 96)
```

This command will give you a list of possible projections for each population, along with number of sites using that projection--format is (projection,sites). Choose either:
* Your original number of individuals (half the total number of possible projections for that population)
* The projection that maximizes the number of sites
  * If doing this, note down this number, since you will this for the 'number of sequences' parameter in your Stairway blueprint file
* It's important to include the -a flag to use all sites; otherwise easySFS will only keep one site per scaffold 

Once you've decided on the appropriate projection(s) (or not), run the following to actually generate the spectr(a/um):
```
python easySFS/easySFS.py -i ./stacks_prelim_test/populations.snps.vcf -p leachs_onepop.popmap -a --proj 224 -o leachs_sfs 
```
Note:
* the --proj flag takes a comma-separated list that equals the number of populations (i.e., if you had three populations, and wanted to project them to 20, 15, and 6 individuals, you'd put --proj 20,15,6. The order of the projection values must match the order of the populations in the --preview; or, you an set this explicitly using the --order flag
* the -o flag sets an output directory for the SFS files
* If you have phased data, you can generate an unfolded SFS. Unfolded spectra are 'better' in that they have more data and can provide higher quality reconstructions. However, in many cases, a folded spectrum (in which you cannot identify derived alleles, so you cut the number of bins in half) is a safer bet. This walkthrough assumes you want to calculate folded spectra
* If you want to generate an unfolded spectrum, just add: --unfolded, and then change the appropriate line in the Stairway blueprint file
* If running on multiple populations, this script will generate both individual population spectrums and also multi-dimensional spectra. The latter can sometimes take a long time, so consider dropping populations before this stage if this ends up being an issue for you.
* If not using the 'vcf-all' flag in populations, use the --total-length flag to specify the total number of genomic sites surveyed, which should help ensure the 0 bin (monomorphic sites) is more accurate

#### Obtain your SFS
Navigate to your output folder, and open the fastsimcoal2 subfolder. View your population of interest's SFS using cat. 
```
cd leachs_sfs/fastsimcoal2
cat LEACH_MAFpop0.obs
```

Sample output
```
1 observation
d0_0    d0_1    d0_2    d0_3    d0_4    d0_5    d0_6    d0_7    d0_8    d0_9    d0_10   d0_11   d0_12   d0_13   d0_14   d0_15   d0_16   d0_17   d0_18   d0_19   d0_20   d0_21d0_22   d0_23   d0_24   d0_25   d0_26   d0_27   d0_28   d0_29   d0_30   d0_31   d0_32   d0_33   d0_34   d0_35   d0_36   d0_37   d0_38   d0_39   d0_40   d0_41   d0_42   d0_43d0_44   d0_45   d0_46   d0_47   d0_48   d0_49   d0_50   d0_51   d0_52   d0_53   d0_54   d0_55   d0_56   d0_57   d0_58   d0_59   d0_60
10556.19919569802 19896.01521123987 10990.62841417877 9359.964578291187 8727.954655170079 8442.922303740483 8293.294725718804 8210.525941897724 8143.49314875933 8087.280338807373 8061.228572287448 8092.004519805983 8263.144439288679 8643.199955645947 9261.586038698679 10054.05570151849 10809.76082850619 11193.56656179125 10949.98432771224 10014.84800341565 8623.941827945506 7091.356912838354 5689.239364009318 4483.74343805637 3625.426216829019 2960.325469944967 2515.424398381165 2321.222349969891 2188.278499242968 2173.368962129355 1081.015098489607 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```
Copy just the numbers (not the d0_0, d0_1 ... etc). Delete the first number (0 bin / monomorphic sites), and omit any bins beyond (1/2 your projection value + 1). You will then paste this set of numbers into your blueprint file under the SFS line below.

To get your L parameter (when using --vcf-all), open R and run...

```
library(tidyverse)
dat = read_table("LEACH_MAFpop0.obs", skip=1)
sum(dat)
```
The number printed should be used as the 'L' parameter in your Stairway blueprint file.

To obtain only the sites you need to copy into easySFS, then run:

```
dat = dat[2:(length(dat)/2 +1)]
print(paste(dat, collapse=" "))
```
Copy this number string (don't copy the quotes), and paste it into the SFS field in your blueprint file.

If you're using an unfolded spectrum, only delete the first bin and last bins (monomorphic sites), then copy all other bins and paste them into your blueprint file under SFS.

### Step 9: Run Stairway Plot v2

Unzip prepared Stairway Plot folder
```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline
unzip stairway_plot_dir.zip
```

I've put a prepared sample blueprint file into this archive. For more information on running Stairway, consult the Stairway GitHub and manuals: https://github.com/xiaoming-liu/stairway-plot-v2/tree/master

Generate blueprint input file
```
cd stairway_plot_v2.1.2
nano seabird_stairway_fold.blueprint
```

Edit this blueprint file to match your species!
Specifically, add...
* popid: The correct popid from your popmap file
* nseq: The number of haploid sequences (i.e., the projection value you chose when running easySFS)
* whether_folded: if using an unfolded spectrum, change this to false
* SFS: take the SFS (without the 0 bin) that you generated in Step 8 and paste it here (should be number of sequences / 2 bins)
* project_dir: the directory you want created to hold your Stairway analysis files and results
* mu: mutation rate--can set this to a dfferent species-specific value if known
* year_per_generation: if you don't know this, can check: https://doi.org/10.1111/cobi.13486
```
#example blueprint file
#input setting
popid: LEACH # id of the population (no white space)
nseq: 224 # number of sequences
L: 10000000 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whether the SFS is folded (true or false)
SFS:    9638.215        3929.77 2243.5499999999997      1493.6750000000002      1110.74 891.8   759.3   667.465 606.0450000000001       567.165 539.245 
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 15 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 7        15      22      28 # number of random break points for each try (separated by white space)
project_dir: seabird_species # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.2e-8 # assumed mutation rate per site per generation
year_per_generation: 24 # assumed generation time (in years)
#plot setting
plot_title: two-epoch_fold # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```
Note:
* Can rename the file if you choose to match each species. If doing so, you'll need to edit the blueprint file name in the following script calls, and also in the SLURM script for running Stairway.

Now, generate input files and script needed to run Stairway plot

```
cd ~/scratch/leachs_storm-petrel_radseq_pipeline/stairway_plot_v2.1.2
module load java
java -cp stairway_plot_es Stairbuilder seabird_stairway_fold.blueprint
```

This command will generate a folder named according to what you put in the project_dir line of the blueprint file. It will also generate a .sh script file that will be named using the name of the blueprint file you used (in this case, seabird_stairway_fold.blueprint.sh)

Finally...actually run Stairway plot! I've included a SLURM script for this in the zipped Stairway archive.
Stairway doesn't usually take a really long time to run, but for species with a large of number of sites, it might take several hours. This script allows 24 hours; can tweak if necessary.

```
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@queensu.ca
#SBATCH --job-name=run_seabird_stairway
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load java/21.0.1

bash seabird_stairway_fold.blueprint.sh
```

```
sbatch run_stairway.sh
```
Note:
* Since this is being run on a compute node, X11 display (which Stairway uses to generate plots) is not loaded, so final plotting will fail. However, the analysis will still run correctly. Once it has finished, you should see a new script file (created by Stairway during its run) which will be called "name_of_your_blueprint_file.blueprint.plot.sh". On a login node (that has X11 enabled), you can now run the final 5 commands from the file.

View the commands using tail
```
tail -n 5 seabird_stairway_fold.blueprint.plot.sh
```
Note:
* If using a larger number of random breakpoints to try (nrand parameter in your blueprint file) adjust the number passed to tail in the -n parameter accordingly

Copy these commands and then run them as below:
```
module load java
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot2 seabird_stairway_fold.blueprint
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot2 seabird_stairway_fold.blueprint rand7
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot2 seabird_stairway_fold.blueprint rand15
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot2 seabird_stairway_fold.blueprint rand22
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot2 seabird_stairway_fold.blueprint rand28
```
This will generate the PDF and PNG plot files.

And you're done!

### Running GONE
* copy GONE folder to your home directory for easiest use later

```
# grab GONE
cd ~
git clone https://github.com/esrud/GONE.git
```

Enter GONE directory, change permissions for all GONE executables

```
cd GONE/Linux/PROGRAMMES/
chmod +x GONE GONEaverage GONEparallel.sh LD_SNP_REAL3 MANAGE_CHROMOSOMES2 SUMM_REP_CHROM3
```
  
### Setup for GONE

* Gone is fairly simple to run, but has a couple idiosyncrasies
* When running Stacks 'populations', need to request plink output files

* To make a 'GONE' format plink file, scaffolds (or chromosomes) need to be renamed to 1-200 (max allowed). Fewer than 200 is fine, especially if you have a higher quality reference genome that's chromosome level. However, plink is human focused, so you'll still need to tweak chromosome naming/etc even if you have a chromosome level assembly, since bird genomes have more chromosomes than humans (especially microchromosomes). Use -99 in the input file to denote that all chromosomes should be used.

* To achieve this, first need to identify the largest 200 scaffolds, retain only those in the plink files, then rename the chromosomes to 1-200 in the map file. If you try to use 'named' chromosomes (i.e., not integers), GONE will crash

```
module load samtools
samtools idxstats <repesentative_bam> > ivory_scaffolds_all.txt 

# sort by length, get top 200 longest scaffolds, write to file
sort -k2 -n -r ivory_scaffolds_all.txt | head -n 200 | cut -f 1 > ivory_scaff_top_200.txt
```

Now, filter plink files to just these largest 200 scaffolds, and do some manual map file editing

```
module load StdEnv/2020 plink/1.9b_6.21-x86_64

# obtain names of all chromosomes to pass into plink
# note: plink requires named chromosomes to be specified --chr <name>,<space><name>, etc...
all_chroms=""
while read p; do all_chroms=`echo -n ${all_chroms}, ${p}`; done < ivory_scaff_top_200.txt
# drop leading comma
all_chroms="${all_chroms:1}"

cd stacks_onepop_3rm
plink --file populations.plink --double-id --allow-extra-chr --chr ${all_chroms} --recode --out ivory_plink_top200_named

cp ivory_plink_top200_named.map ivory_plink_top200.map
cp ivory_plink_top200_named.ped ivory_plink_top200.ped

# get new chrom list in case any chromosomes had no variants
cat ivory_plink_top200.map | cut -f 1 | sort -u -k 1 > new_chrom_list_200.txt

# make look-up file
touch chrom_200_lookup.txt
var_i=1
while read p; do echo ${p},${var_i} >> chrom_200_lookup.txt; var_i=`expr $var_i + 1`; done < new_chrom_list_200.txt

# convert chromosome names to numbers in map file (for GONE formatting)
while read p; do line1=`echo $p | cut -f 1 -d ","`; line2=`echo $p | cut -f 2 -d ","`; sed -i "s/$line1/$line2/g" ivory_plink_top200.map; done < chrom_200_lookup.txt
```

Now, create a working directory for your GONE analysis

```
# create a working directory for your GONE analysis
mkdir ivgu_gone
cd ivgu_gone

# copy in GONE executables and scripts
cp -r ~/GONE/Linux/INPUT_PARAMETERS_FILE ~/GONE/Linux/PROGRAMMES/ ~/GONE/Linux/script_GONE.sh .

# copy in your map and ped files
cp ~/<your_working_dir>/data.map ~/<your_working_dir>/data.ped
```

Now, do any necessary editing of INPUT_PARAMETER_FILE


Run GONE
* actually runs quite quickly, can run in a salloc

```
salloc --mem 40G --time 2:00:00

bash script_GONE.sh <file_prefix_before.ped/.map>
```

Looping GONE over multiple colonies

* You will need your colony popmap file (my file was popmap_adults.txt)

This first (long) loop will generate one GONE format plink file for every unique population in your popmap. You will need to tweak paths to your plink files to match your directories, and will need to tweak the 'family' portion of the keep file.

Other notes...
* Run this after you've already run GONE initially (makes use of the scaffold files, and the 'named' top200 plink files). If you've already used my script to convert your top200 into GONE format, you'll unfortunately have to re-run the preceding step
* You may need to check your 'family' names (actually your popmap names). The files I generated are for my specific ivory gull situation, where I had named the overall population 'global'. Whatever you named your population should go in place of global (for whichever population you're splitting)
* change the paths of the plink files to match your directory tree

```
for POP in `cut -f 2 popmap_adults.txt | sort -u -k 1`; do grep $POP popmap_adults.txt | cut -f 1 > split_popmap_${POP}.txt; for LINE in `cat split_popmap_${POP}.txt`; do touch temp.txt; printf "${POP}\t${LINE}\n" >> temp.txt; done; mv temp.txt split_popmap_${POP}.txt; plink --file stacks_onepop_3rm/ivory_plink_top200_named --keep split_popmap_${POP}.txt --double-id --allow-extra-chr --recode --out stacks_onepop_3rm/ivory_plink_top200_${POP}; while read p; do line1=`echo $p | cut -f 1 -d ","`; line2=`echo $p | cut -f 2 -d ","`; sed -i "s/$line1/$line2/g" stacks_onepop_3rm/ivory_plink_top200_${POP}.map; done < stacks_onepop_3rm/chrom_200_lookup.txt; done

```

Open the folder containing your plink files, and make a new file called "run_gone.sh". The script file should contain:

```
# create a working directory for your GONE analysis
wd=$1
map=$2
ped=$3
prefix=$4

mkdir $wd
cd $wd

# copy in GONE executables and scripts
cp -r ~/GONE/Linux/INPUT_PARAMETERS_FILE ~/GONE/Linux/PROGRAMMES/ ~/GONE/Linux/script_GONE.sh .

# copy in your map and ped files
cp $map $ped .

bash script_GONE.sh $prefix
```

Make the file executable

```
chmod +x run_gone.sh
```

Now, open a salloc and run GONE once for every plink POP file
* again, will need to tweak paths. Make sure that the paths to your map and ped files are absolute paths!
* If you have a lot of populations, may need to increase the amount of time requested
```
salloc --mem 40G --time 2:00:00

# change the absolute paths here to the appropriate ones for your cluster (to the folder containing your plink files)
for POP in `cut -f 2 ../popmap_adults.txt | sort -u -k 1`; do bash run_gone.sh gone_${POP} /global/home/hpc5400/scratch/ivory_gull_stairway/leachs_storm-petrel_radseq_pipeline/stacks_onepop_3rm/ivory_plink_top200_${POP}.map /global/home/hpc5400/scratch/ivory_gull_stairway/leachs_storm-petrel_radseq_pipeline/stacks_onepop_3rm/ivory_plink_top200_${POP}.ped ivory_plink_top200_${POP}; done

```

Preliminary R plotting code

```
library(tidyverse)
dat = read_delim("Output_Ne_ivory_plink_top200", skip=1)

# plot in units of generation time
pdf("base_r_gone_testplot_ivgu.pdf")
plot(Geometric_mean~Generation, data=dat)
dev.off()

# convert to 'actual' time
dat$raw_time = dat$Generation * 7.9
pdf("base_r_gone_testplot_ivgu_time.pdf", width=18, height=8)
plot(Geometric_mean~raw_time, xlab="Time",  data=dat)
dev.off()
```

