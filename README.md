## Leach's storm-petrel ddRadSeq pipeline
Prepared by Chris Boccia, Nov. 2024


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

populations -P ./stacks_both -O ./stacks_prelim_test -r 0.75 -R 0.5 --min-maf 0.05 -M popmap_w_batches.txt -t 8 --radpainter --vcf --genepop --structure --plink --phylip --treemix
```

```
sbatch populations_run.sh
```
Notes:
* This script will generate input formats for a wide variety of downstream programs. VCF files are the most widely used, but radpainter, structure, plink, and genepop are also commonly used. Phylip and treemix are mainly used for explicit phylogenetic reconstruction and comparative methods.

### Steps before downstream processing

At this point, most folks use easySFS to generate a site frequency spectrum, and plink or similar to do PCA

Notes on using Stairway plot:
* To correctly estimate Ne, you need to provide Stairway plot with an accurate 'L'. The L parameter is the length of sequence/number of sites surveyed for variation. This changes depending on your filtering, and Ne can be over/under estimated if you get this parameter wrong. To determine what it should be for your data set, following your final filtering using populations, run:

```
tail -n  20  populations.log
```

and look for the 'genomic sites' parameter reported here:

```
Removed 1060001 loci that did not pass sample/population constraints from 1102457 loci.
Kept 42456 loci, composed of 13684652 sites; 6246148 of those sites were filtered, 62363 variant sites remained.
    7577896 genomic sites, of which 14449 were covered by multiple loci (0.2%).
Mean genotyped sites per locus: 178.85bp (stderr 0.42).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  Baccalieu: 34.057 samples per locus; pi: 0.29648; all/variant/polymorphic sites: 7493660/61235/60938; private alleles: 0
  Gull: 43.048 samples per locus; pi: 0.29775; all/variant/polymorphic sites: 7481718/60974/60867; private alleles: 0
  Country: 22.265 samples per locus; pi: 0.29624; all/variant/polymorphic sites: 7402803/60971/59778; private alleles: 0
  Kent: 21.227 samples per locus; pi: 0.29657; all/variant/polymorphic sites: 7446907/60919/59605; private alleles: 0
  Bon_Portage: 19.162 samples per locus; pi: 0.29855; all/variant/polymorphic sites: 7457632/60866/59410; private alleles: 0
  Middle_Lawn: 16.749 samples per locus; pi: 0.29195; all/variant/polymorphic sites: 7454154/58622/56357; private alleles: 0
  Corossol: 15.374 samples per locus; pi: 0.29624; all/variant/polymorphic sites: 7461236/61066/58292; private alleles: 0
  Hernyken: 23.063 samples per locus; pi: 0.29393; all/variant/polymorphic sites: 7566670/60947/59936; private alleles: 0
  Unknown: 79.462 samples per locus; pi: 0.29843; all/variant/polymorphic sites: 7564697/61886/61879; private alleles: 0
  Green: 20.407 samples per locus; pi: 0.29744; all/variant/polymorphic sites: 7474731/61594/60324; private alleles: 0
  Baja: 3 samples per locus; pi: 0.14971; all/variant/polymorphic sites: 7355412/55459/18050; private alleles: 0
  Iceland: 14.681 samples per locus; pi: 0.29676; all/variant/polymorphic sites: 7535137/60468/58054; private alleles: 0
Populations is done.
```

I believe this *should* be the correct L to use, as it includes both monomorphic and polymorphic sites.

<br>

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
python easySFS/easySFS.py -i ./stacks_prelim_test/populations.snps.vcf -p leachs_onepop.popmap --preview
```

Sample output (from Leach's data):
```
LEACH
(2, 117)        (3, 175)        (4, 212)        (5, 240)        (6, 261)        (7, 278)        (8, 293)        (9, 305)        (10, 315)       (11, 324)       (12, 331)   (13, 338)        (14, 344)       (15, 350)       (16, 354)       (17, 359)       (18, 362)       (19, 366)       (20, 369)       (21, 372)       (22, 375)       (23, 377)   (24, 379)        (25, 381)       (26, 383)       (27, 385)       (28, 386)       (29, 387)       (30, 389)       (31, 390)       (32, 391)       (33, 392)       (34, 393)   (35, 394)        (36, 395)       (37, 395)       (38, 396)       (39, 397)       (40, 397)       (41, 398)       (42, 398)       (43, 399)       (44, 399)       (45, 399)   (46, 400)        (47, 400)       (48, 401)       (49, 401)       (50, 401)       (51, 401)       (52, 402)       (53, 402)       (54, 402)       (55, 402)       (56, 402)   (57, 403)        (58, 403)       (59, 403)       (60, 403)       (61, 403)       (62, 403)       (63, 403)       (64, 403)       (65, 404)       (66, 404)       (67, 404)   (68, 404)        (69, 404)       (70, 404)       (71, 404)       (72, 404)       (73, 404)       (74, 404)       (75, 404)       (76, 404)       (77, 404)       (78, 404)   (79, 404)        (80, 404)       (81, 405)       (82, 405)       (83, 405)       (84, 405)       (85, 405)       (86, 405)       (87, 405)       (88, 405)       (89, 405)   (90, 405)        (91, 405)       (92, 405)       (93, 405)       (94, 405)       (95, 405)       (96, 405)       (97, 405)       (98, 405)       (99, 405)       (100, 405)  (101, 405)       (102, 405)      (103, 405)      (104, 405)      (105, 405)      (106, 405)      (107, 405)      (108, 405)      (109, 405)      (110, 405)      (111, 405)  (112, 405)       (113, 405)      (114, 405)      (115, 405)      (116, 405)      (117, 405)      (118, 405)      (119, 405)      (120, 405)      (121, 405)      (122, 405)  (123, 405)       (124, 405)      (125, 405)      (126, 405)      (127, 405)      (128, 405)      (129, 405)      (130, 405)      (131, 405)      (132, 405)      (133, 405)  (134, 405)       (135, 405)      (136, 405)      (137, 405)      (138, 405)      (139, 405)      (140, 405)      (141, 405)      (142, 405)      (143, 405)      (144, 405)  (145, 405)       (146, 405)      (147, 405)      (148, 405)      (149, 405)      (150, 405)      (151, 405)      (152, 405)      (153, 405)      (154, 405)      (155, 405)  (156, 405)       (157, 405)      (158, 405)      (159, 405)      (160, 405)      (161, 405)      (162, 405)      (163, 405)      (164, 405)      (165, 405)      (166, 405)  (167, 405)       (168, 405)      (169, 405)      (170, 405)      (171, 405)      (172, 405)      (173, 405)      (174, 405)      (175, 405)      (176, 405)      (177, 404)  (178, 404)       (179, 404)      (180, 404)      (181, 404)      (182, 404)      (183, 404)      (184, 404)      (185, 402)      (186, 402)      (187, 402)      (188, 402)  (189, 402)       (190, 402)      (191, 402)      (192, 402)      (193, 401)      (194, 401)      (195, 401)      (196, 401)      (197, 401)      (198, 401)      (199, 401)  (200, 401)       (201, 400)      (202, 400)      (203, 400)      (204, 400)      (205, 400)      (206, 400)      (207, 399)      (208, 399)      (209, 399)      (210, 399)  (211, 397)       (212, 397)      (213, 396)      (214, 396)      (215, 396)      (216, 396)      (217, 396)      (218, 396)      (219, 396)      (220, 396)      (221, 395)  (222, 395)       (223, 395)      (224, 395)      (225, 395)      (226, 395)      (227, 394)      (228, 394)      (229, 394)      (230, 394)      (231, 394)      (232, 394)  (233, 394)       (234, 394)      (235, 394)      (236, 394)      (237, 394)      (238, 394)      (239, 394)      (240, 394)      (241, 394)      (242, 394)      (243, 393)  (244, 393)       (245, 391)      (246, 391)      (247, 390)      (248, 390)      (249, 390)      (250, 390)      (251, 390)      (252, 390)      (253, 390)      (254, 390)  (255, 389)       (256, 389)      (257, 387)      (258, 387)      (259, 387)      (260, 387)      (261, 387)      (262, 387)      (263, 386)      (264, 386)      (265, 386)  (266, 386)       (267, 386)      (268, 386)      (269, 386)      (270, 386)      (271, 386)      (272, 386)      (273, 385)      (274, 385)      (275, 385)      (276, 385)  (277, 384)       (278, 384)      (279, 384)      (280, 384)      (281, 384)      (282, 384)      (283, 384)      (284, 384)      (285, 384)      (286, 384)      (287, 382)  (288, 382)       (289, 382)      (290, 382)      (291, 379)      (292, 379)      (293, 379)      (294, 379)      (295, 379)      (296, 379)      (297, 379)      (298, 379)  (299, 379)       (300, 379)      (301, 379)      (302, 379)      (303, 379)      (304, 379)      (305, 377)      (306, 377)      (307, 375)      (308, 375)      (309, 375)  (310, 375)       (311, 374)      (312, 374)      (313, 374)      (314, 374)      (315, 374)      (316, 374)      (317, 374)      (318, 374)      (319, 373)      (320, 373)  (321, 373)       (322, 373)      (323, 371)      (324, 371)      (325, 369)      (326, 369)      (327, 368)      (328, 368)      (329, 367)      (330, 367)      (331, 366)  (332, 366)       (333, 365)      (334, 365)      (335, 364)      (336, 364)      (337, 363)      (338, 363)      (339, 363)      (340, 363)      (341, 361)      (342, 361)  (343, 359)       (344, 359)      (345, 359)      (346, 359)      (347, 356)      (348, 356)      (349, 355)      (350, 355)      (351, 353)      (352, 353)      (353, 350)  (354, 350)       (355, 350)      (356, 350)      (357, 349)      (358, 349)      (359, 347)      (360, 347)      (361, 346)      (362, 346)      (363, 345)      (364, 345)  (365, 343)       (366, 343)      (367, 341)      (368, 341)      (369, 338)      (370, 338)      (371, 336)      (372, 336)      (373, 336)      (374, 336)      (375, 335)  (376, 335)       (377, 334)      (378, 334)      (379, 333)      (380, 333)      (381, 332)      (382, 332)      (383, 330)      (384, 330)      (385, 329)      (386, 329)  (387, 328)       (388, 328)      (389, 328)      (390, 328)      (391, 327)      (392, 327)      (393, 325)      (394, 325)      (395, 322)      (396, 322)      (397, 321)  (398, 321)       (399, 318)      (400, 318)      (401, 316)      (402, 316)      (403, 313)      (404, 313)      (405, 310)      (406, 310)      (407, 305)      (408, 305)  (409, 302)       (410, 302)      (411, 296)      (412, 296)      (413, 294)      (414, 294)      (415, 289)      (416, 289)      (417, 286)      (418, 286)      (419, 282)  (420, 282)       (421, 268)      (422, 268)      (423, 258)      (424, 258)      (425, 250)      (426, 250)      (427, 234)      (428, 234)      (429, 219)      (430, 219)  (431, 203)       (432, 203)      (433, 185)      (434, 185)      (435, 162)      (436, 162)      (437, 135)      (438, 135)      (439, 110)      (440, 110)      (441, 75)   (442, 75)        (443, 52)       (444, 52)       (445, 22)       (446, 22)       (447, 4)        (448, 4)
```

This command will give you a list of possible projections for each population, along with number of sites using that projection--format is (projection,sites). Choose either:
* Your original number of individuals (half the total number of possible projections for that population)
* The projection that maximizes the number of sites
  * If doing this, note down this number, since you will this for the 'number of sequences' parameter in your Stairway blueprint file

Once you've decided on the appropriate projection(s) (or not), run the following to actually generate the spectr(a/um):
```
python easySFS/easySFS.py -i ./stacks_prelim_test/populations.snps.vcf -p leachs_onepop.popmap --proj 224 -o leachs_sfs
```
Note:
* the --proj flag takes a comma-separated list that equals the number of populations (i.e., if you had three populations, and wanted to project them to 20, 15, and 6 individuals, you'd put --proj 20,15,6. The order of the projection values must match the order of the populations in the --preview; or, you an set this explicitly using the --order flag
* the -o flag sets an output directory for the SFS files
* If you have phased data, you can generate an unfolded SFS. Unfolded spectra are 'better' in that they have more data and can provide higher quality reconstructions. However, in many cases, a folded spectrum (in which you cannot identify derived alleles, so you cut the number of bins in half) is a safer bet. This walkthrough assumes you want to calculate folded spectra
* If you want to generate an unfolded spectrum, just add: --unfolded, and then change the appropriate line in the Stairway blueprint file
* If running on multiple populations, this script will generate both individual population spectrums and also multi-dimensional spectra. The latter can sometimes take a long time, so consider dropping populations before this stage if this ends up being an issue for you.

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
Copy just the numbers (not the d0_0, d0_1 ... etc). Delete the first number (0 bin). You will paste this into your blueprint file under the SFS line below.

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
* nseq: The number of sequences (number of individuals * 2, since they're diploid)
* L (number of observed sites, including monomorphic sites). See text in previous sections.
* whether_folded: if using an unfolded spectrum, change this to false
* SFS: take the SFS (without the 0 bin) that you generated in Step 8 and paste it here
* project_dir: the directory you want created to hold your Stairway analysis files and results
* mu: mutation rate--can set this to a dfferent species-specific value if known
* year_per_generation: if you don't know this, can check: https://doi.org/10.1111/cobi.13486
```
#example blueprint file
#input setting
popid: LEACH # id of the population (no white space)
nseq: 30 # number of sequences
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
* Since this is being run on a compute node, X11 display (which Stairway uses to generate plots) is not loaded, so final plotting will fail. However, the analysis will still run correctly. Once it has finished, you should see a new script file (created by Stairway during its run) which will be called "name_of_your_blueprint_file.blueprint.plot.sh". On a login node (with a salloc), you can now run the final 5 commands from the file.

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




