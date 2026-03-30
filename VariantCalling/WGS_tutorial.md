# Whole genome sequencing analysis tutorial
In this tutorial we will go over how to analyse Illumina WGS data. This will comprise four main steps:

1. Sequencing data QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read trimming ([fastp](https://github.com/FelixKrueger/TrimGalore/))
3. Read alignment with substeps ([BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2))
   1. Indexing reference
   2. Mapping the reads
   3. Sorting, compressing and indexing
   3. Gather metrics
4. Mark duplicates ([fastdup]())
5. Call variants ([FreeBayes](https://github.com/freebayes/freebayes))
6. Annotate variants ([VEP](https://www.ensembl.org/info/docs/tools/vep/))
7. Collect all metrics ([MultiQC](https://multiqc.info/docs/))

### Install dependencies
Before starting the analyses, we need to install all the tools required.
To do so, we will use [miniforge](), that enables us to install several dependencies easily, without the need to compile anything:
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```
And then run it with:
```
bash Miniforge3-$(uname)-$(uname -m).sh
```
We follow the guided installation to completely install miniforge. Then, we install all the dependencies with the command:
```
mamba create -n readmapping_env -c conda-forge -c bioconda samtools multiqc vcftools bcftools bwa-mem2 fastdup fastp mosdepth fastqc gatk
mamba create -n variant_calling_env -c conda-forge -c bioconda clair3 clair3-illumina
mamba create -n annotation_env -c conda-forge -c bioconda ensembl-vep 
```
We need to create multiple environments to avoid issues with the dependencies.


### Resources
Before proceeding with the analyses, we need to download the raw sequencing data. The files are in the same format described for the RNA-seq analyses, and comprise of two compressed FASTQ files.

First, we create the folder that will contain the data using `mkdir`:
```
mkdir DATA && cd DATA
```

We can now download the sequencing data into the new folder. To do this we will use the [wget](https://en.wikipedia.org/wiki/Wget) program that can pull files from the internet (in this case dropbox):

```
wget https://www.dropbox.com/s/2nw3bokejpueudp/Holstein_1_R1.fq.gz
wget https://www.dropbox.com/s/zgmxmxkwvtkk5qf/Holstein_1_R2.fq.gz
```

Note there are two sets of reads as this is [paired end](https://emea.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) data.

Finally, we exit from the `DATA` folder:
```
cd ../
```

## Preprocessing the reads
### Quality check of the reads
Similarly to what done for the RNA-seq data, the first thing to do in a sequencing project is to check the quality of the data.
We can check the quality of Illumina data is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) it for all our samples
using a `for` loop.
First, we create the output directory named `FASTQC` with `mkdir`:
```
mkdir FASTQC
```

Then, for each sample name present in the DATA folder, we can 
```
fastqc DATA/Holstein_1_R1.fq.gz DATA/Holstein_1_R2.fq.gz -o FASTQC/ --threads 1 
```

This will generate reports in HTML (browser format) that describe the quality of the reads.

Have a look at your FastQC reports.

> *How many sequence reads are there in each file?*

> *Which metrics do you have warnings on?*


### Trimming the reads
Following the computation of the FastQC reports, we can trim our reads to remove all low-quality bases and sequences.
We will run `trim_galore`, an popular utility able to automatically identify and clean our raw reads.
First, we create the output folder:
```
mkdir TRIM
```
Then we run trim_galore on each fastq pair:
```
fastp \
   --in1 DATA/Holstein_1_R1.fq.gz \
   --in2 DATA/Holstein_1_R2.fq.gz \
   --out1 TRIM/Holstein_1_R1.trimmed.fq.gz \
   --out2 TRIM/Holstein_1_R2.trimmed.fq.gz \
   --html TRIM/Holstein_1.fastp.html \
   --json TRIM/Holstein_1.fastp.json \
   --thread 4 \
   --detect_adapter_for_pe
```
This will generate a set of trimmed fastq files for each sample, as well as generate the new FastQC 
reports for the output datasets. 

> *How many bases have been trimmed for quality in the R1 library?*


## Alignment of the reads
### Prepare the reference genome
Reads alignment involves a series of step 
Before running the alignment, we need to download and prepare the reference genome to process.
We first create the folder where to save our reference genome:
```
mkdir REFERENCE && cd REFERENCE
```

To facilitate the process, we provide a reference sequence comprising of only chromosome 20:
```
wget https://www.dropbox.com/s/dq2qgdko6pgk14o/chr20.fa
```
And then we can exit the folder:
```
cd ..
```

Firstly, we prepare the `fai` index, which summarises the sequence names, size and position in the file.
We can generate this using `samtools`, a popular tool used to process sequence data, with the `faidx` command:
```
samtools faidx REFERENCE/chr20.fa
```
This will generate a file called `REFERENCE/chr20.fa.fai`, which contains the information concerning the sequence.
In particular, the name of the sequences will be in the first column, whereas the size of the sequence will be in the second column.

> *How many bases are present in chromosome 20?*

Then, we prepare the indexes needed from BWA to align the reads of interest:
```
bwa-mem2 index REFERENCE/chr20.fa
```

This process will take seconds to minutes, and will generate a set of large binary files. 
These files are not meant to be read and used by humans, but will be used by BWA to align
the reads efficiently.

### Align the reads
The next step required to call variants is to align the reads to the reference genome.
To do so we will use bwa-mem2, a faster and more efficient version of the standard bwa-mem 
aligner.
We first create the directory where to save out results:
```
mkdir ALIGN
```

Then, we align each sample and save it in a separate subfolder
```
bwa-mem2 mem -R "@RG\tID:Holstein_1\tSM:Holstein_1" -t 1 REFERENCE/chr20.fa TRIM/Holstein_1_R1_val_1.fq.gz TRIM/Holstein_1_R2_val_2.fq.gz > ALIGN/Holstein_1.raw.sam
```

This will save the aligned reads in SAM (Sequence Alignment Map) format.


### Sorting and indexing the alignments
Following the alignments, the reads will be collected in an unordinate way in the output files.
Unsorted reads represent a relevant bottleneck in the downstream analyses, since accessing ordinate files
is quicker.
To sort the files, we are going to use [samtools]() one more:
```
samtools sort -O BAM -@1 ALIGN/Holstein_1.raw.sam > ALIGN/Holstein_1.sort.bam
```

With this command, we are doing two separate operations:
 1. Sort the input sam file
 2. Compress the sam to bam (binary alignment mapping) to save space

Following the sorting, we can index the file. Indexing is an operation that allows to other tools to 
access specific portions of a file without loading it into the memory as a whole or going through it 
looking for the desired region.
To index, we resort to samtools once more:
```
samtools index ALIGN/Holstein_1.sort.bam
```

This will generate files ending with the `.bai` suffix, which are the indexes of each bam file.

### Compute statistics on the alignments
After aligning the reads to the reference, we can gather some metrics on the alignments.
To do so, we use samtools again to compute both alignments metrics (`samtools flagstat`) and coverage (`samtools coverage`) of the genome:
```
samtools flagstat ALIGN/Holstein_1.sort.bam > ALIGN/Holstein_1.flagstat
samtools coverage ALIGN/Holstein_1.sort.bam > ALIGN/Holstein_1.coverage
```
The `flagstat` file contains the key metrics for the alignments, including number of aligned reads,
primary and secondary alignments and more. The second collects the details of the coverage of the sam file.

> *How many properly paired reads are present in the file?*
> *How many secondarily aligned reads are present in the file?*
> *How many bases are covered (i.e. have at least one alignment) in the coverage file?*

## Alignments post-processing
### Filtering the reads
Following the alignment, we can filter our reads for a number of parameters. There are a number of 
parameters that can be used to filter out the reads, and can be seen [here](https://broadinstitute.github.io/picard/explain-flags.html).
First, we create a new folder for the filtered reads:
```
mkdir FILTER
```

Then we perform the actual filtering. In this example, we are going to drop all read with mapping quality 
equal to 0 and unmapped reads:
```
samtools view -h -b -F 4 -q 1 ALIGN/Holstein_1.sort.bam > FILTER/Holstein_1.filt.bam
samtools index FILTER/Holstein_1.filt.bam
```
The option provided are the following:
 1. `-F 4`: filter out all unmapped reads
 2. `-q 1`: drop reads with mapping quality below 1
 3. `-h`  : save the header of the bam file in the output
 4. `-b`  : save the output as bam file


### Remove duplicates
The raw alignment files may contain artificial PCR duplicates due to the chemistry required to prepare the samples.
Being identical, these reads might ...

First, we create a new folder where to save our deduplicated files:
```
mkdir MARKDUP
```

Then, we mark duplicated reads by using `fastdup` as follow:
```
fastdup 
```
This command will run `fastdup` on the input bam file, saving an output bam file in the 
output folder specified using 4 threads for job.

> *What is the percentage of duplicated reads?*

## Variant calling
Following the generation of the mark duplicates bam file, we can progress to perform the actuall 
variant calling. This process can be performed using a number of possible software solutions.
Among the most popular, we have to mention:
1. [bcftools](https://samtools.github.io/bcftools/bcftools.html) mpileup + call: This is the most classical of the variant calling method. Good with SNPs, less so with the InDels
2. [GATK](gatk.broadinstitute.org/) HaplotypeCaller: Golden standard for a long time, it performs a local re-assembly of the reads to better resolve InDels
3. [FreeBayes](https://github.com/freebayes/freebayes): another commonly adopted variant caller, uses bayesian inference to in place of the maximum likelyhood method from 
4. [DeepVariant](https://github.com/google/deepvariant): Google's own variant caller, uses an AI image classification approach to quickly identify the variants
5. [Clair3](https://github.com/HKU-BAL/Clair3): another toold using AI and Deep Learning to detect variants in the reads provided

In our case, we are going to use Clair3, which is a popular, well-established variant caller with plenty of 
functionality. Clair3 has a large number of parameters that can be fine-tuned to achieve greater sensitivity, 
to improve the calling, or to processed pooled samples at once. 

To run Clair3 we need to first download the appropriate model for the illumina reads. The officially trained models can be found [here](https://github.com/HKU-BAL/Clair3#HKU-provided-Models). We need to download the illumina one as follow:
```
mkdir -p models/ilmn
wget -O models/ilmn/pileup.pt https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/ilmn/pileup.pt
wget -O models/ilmn/full_alignment.pt https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/ilmn/full_alignment.pt
```

To perform the variant calling we are going to create a new folder where to save our output VCFs:
```
mkdir CLAIR3
```

Then, we can proceed with the variant calling. In our case, we are going to run the caller with base parameters.
We can call the variants as follow: 
```
run_clair3.sh \
  --bam_fn=MARKDUP/Holstein_1.md.bam \ ## change your bam file name here
  --ref_fn=REFERENCE/chr20.fa \        ## change your reference file name here
  --threads=4 \                        ## maximum threads to be used
  --platform="ilmn" \                  ## options: {ont,hifi,ilmn}
  --model_path="${PWD}/models/ilmn" \  ## Path to the model
  --gvcf \                             ## Create also the GVCF file
  --output=CLAIR3/                     ## output directory
```
This will generate an output VCF file ready for downstream analyses.

## Compressing and indexing the VCF files
It is best practice to compress and index the VCF files generated. The compression allows to save space,
while the indexing speed up the access to the entries of the file. We will compress the VCF file using 
[bgzip](https://github.com/samtools/tabix), a compression algorithm optimized to work with VCF files:
```
bgzip FBAYES/Holstein_1.vcf
```
This will replace our VCF file with a `.vcf.gz`, a compressed version of this.

Then, we will index the file using [tabix](https://github.com/samtools/tabix):
```
tabix -p vcf FBAYES/Holstein_1.vcf.gz
```
This will create a file in `.tbi` format, in which the index will be stored.

An alternative method to index is by using the `bcftools index` command:
```
bcftools index -t FBAYES/Holstein_1.vcf.gz
```

## VCF metrics and statistics
The VCF files generated by any variant caller are large files that store a large amount of information.
This makes them very hard to screen by hand, and not easy to visualise. To have a peek to the file themself, we 
can use the `less` command:
```
less FBAYES/Holstein_1.vcf.gz
```
Since the lines of the file are long, it can be convenient to add the `-S` option to `less`, which avoid the automatic new line for each row
in the VCF:
```
less -S FBAYES/Holstein_1.vcf.gz
```
In this case, we can use all the arrows to move around the file vertically (up/down) and horizontally (right/left).

Given the large size of these files, we can use some software to collect metrics to describe the file.
One of the commonly used tools is [bcftools](https://samtools.github.io/bcftools/bcftools.html), which can screen the VCF for us
and print out generic metrics: 
```
bcftools stats FBAYES/Holstein_1.vcf.gz > FBAYES/Holstein_1.stats
```
Remember, the `>` character means "**redirection**": it redirect the output of a command from the screen (known as STDOUT) to a 
file of destination.


## Filter variants
Combining the files generates a VCF file with a large nubmer fo entries. Several of these might have low coverage
issues, low quality or high missingness rate.
These can be filtered out using [vcftools](vcftools.github.io/). We first exclude variants with low QUAL score, 
suggesting that these are likely wrong calls:
```
vcftools --gzvcf FBAYES/Holstein_1.vcf.gz --minQ 20 --recode --recode-INFO-all --out FBAYES/Holstein_1.highQ
```

Then, we filter out variants with too low depth of sequencing. The depth of sequencing is the number of reads aligned to a given base,
and carrying either the reference or alternative allele. Regions with too low depth of sequencing are more difficult to distinguish
from errors, whereas regions with too high depth are likely to be repetitive regions, enriched with sequencing errors.
The VCF file codify the depth of sequencing in the `DP` field, and is associated to each genotype for each individual.  
To filter out variants with low (depth < 5) or high (depth > 50) depth of sequencing, we use `vcftools` again:
```
vcftools --vcf FBAYES/Holstein_1.highQ.recode.vcf --min-meanDP 5 --max-meanDP 60 --recode --recode-INFO-all --out FBAYES/Holstein_1.highQ.5-60DP
```

Additional possible filtering are accessible in the [vcftools documentation](https://vcftools.github.io/man_latest.html).

> *How many variants do we save after the filtering are applied?*

## Variant Effect Predictor
Following the quick filtering of the VCF file, we can annotate the variants using the Variant Effect Predictor ([VEP](https://www.ensembl.org/info/docs/tools/vep)) tool.
VEP is a complex software, with multiple options available and configuration that can be used. It also undergoes to frequent revision to incorporate new data, making
it important to point out which version of the software we are using. In our case, we are processing the data using the latest version of the software (v106.1 at the 
time of writing this tutorial). 
First off, we create the output directory:
```
mkdir VEP
```

Then, we annotate the vcf file using the VEP software:
```
vep -i FBAYES/Holstein_1.highQ.5-60DP.recode.vcf --species bos_taurus -o VEP/Holstein_1.highQ.5-60DP.recode.vep.vcf --vcf --database
```

This command generate two separate outputs:
1. An annotated VCF file
2. An HTML file with the report of the annotated variants

> *How many variants with moderate effect do we have?*

## Create a summary report with MultiQC
Finally, we can collect generic metrics for the whole analysis using [MultiQC](https://multiqc.info/docs/), sinmilarly to what done for the RNA-seq analyses.
To do so, we can simply run:
```
multiqc .
```
The resulting report will be saved as `multiqc_report.html`, and can be visualized using any browser.