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

* If using on a different dataset, you need to sub in your fastq file names (R1 and R2) into each unparallelized DBR scrip (dbr01_trim.sh, dbr08_trim.sh, dbr10_trim.sh, drb11_trim.sh)

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
* This script will launch a separate job for each of the 8 different paired chunks of the batch 2 fastq file (i.e., x##_R1.fastq and X##_R2.fastq)
* The --array-10-17 slurm flag causes the SLURM_ARRAY_TASK_ID variable to take on the numbers between 10 and 17, which match our split files (x10-17)
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
sbatch align_bwa.sh
```

Notes:
* As with the trimming script, this script will output alignment files to the same directory as the source trimmed fastq files
* This step combines forward and reverse reads, so you get only a single output file in bam format
* The script also sorts each file, which I believe is necessary for Stacks (or at least makes running Stacks faster)

First, since some samples overlap between batches, we need to add a prefix to each bam file

Batch 1
```
cd batch1/batch1_demultiplexed

# for each bam file, add b1_ to the start of the file name
for FILE in `ls -1 *.bam`; do
mv $FILE b1_${FILE}
done
```

Batch 2
```
cd batch2/batch2_demultiplexed

# for each bam file, add b1_ to the start of the file name
for FILE in `ls -1 *.bam`; do
mv $FILE b1_${FILE}
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

Run the following script to generate a Stacks output folder. You can then use this output folder to generate filtered datasets; filtering is very quick and lightweight, so you can do it without using a SLURM script.

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

### Downstream processing

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
