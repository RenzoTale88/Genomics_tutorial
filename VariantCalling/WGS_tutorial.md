# Whole genome sequencing analysis tutorial
In this tutorial we will go over how to analyse Illumina WGS data. This will comprise four main steps:

1. Sequencing data QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
1. Read trimming ([fastp](https://github.com/FelixKrueger/TrimGalore/))
1. Read alignment with substeps ([BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2))
   1. Indexing reference
   1. Mapping the reads
   1. Sorting, compressing and indexing
   1. Gather metrics
1. Mark duplicates ([fastdup](https://github.com/zzhofict/FastDup))
1. Call variants ([Clair3](https://github.com/HKU-BAL/Clair3))
1. Combine the G.VCF ([GLNexus](https://github.com/dnanexus-rnd/GLnexus))
1. Annotate variants ([snpEff](https://pcingola.github.io/SnpEff/))
1. Collect all metrics ([MultiQC](https://multiqc.info/docs/))

### Install dependencies
Before starting the analyses, we need to install all the tools required.
For this tutorial we will use [anaconda](https://anaconda.org/),an open source distribution platform initially designed for python packages, but that can ship also other tools. Anaconda works with [virtual environments](https://en.wikipedia.org/wiki/Virtual_environment_software), a contained isolated installations and dependencies.
To do so, we will use [miniforge](https://github.com/conda-forge/miniforge), that enables us to install several dependencies easily, without the need to compile anything:
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```
And then run it with:
```
bash Miniforge3-$(uname)-$(uname -m).sh
```

We then run:
```
source ~/.bashrc
```

We follow the guided installation to completely install miniforge. Then, we install all the dependencies with the command:
```
mamba create -n readmapping_env -c conda-forge -c bioconda samtools multiqc vcftools bcftools bwa-mem2 fastdup fastp mosdepth fastqc gatk
mamba create -n variant_calling_env -c conda-forge -c bioconda clair3 clair3-illumina glnexus bcftools vcftools snpeff
```
We need to create multiple environments to avoid issues with the dependencies.


### Resources
Before proceeding with the analyses, we need to download the raw sequencing data. The files are in the same format described for the RNA-seq analyses, and comprise of two compressed FASTQ files.

First, we create the folder that will contain the data using `mkdir`:
```
mkdir DATA
```

We can now download the sequencing data into the new folder. To do this we will use the [wget](https://en.wikipedia.org/wiki/Wget) program that can pull files from the internet (in this case dropbox):

```
wget -O DATA/Holstein_1_R1.fq.gz https://www.dropbox.com/s/2nw3bokejpueudp/Holstein_1_R1.fq.gz
wget -O DATA/Holstein_1_R2.fq.gz https://www.dropbox.com/s/zgmxmxkwvtkk5qf/Holstein_1_R2.fq.gz
```

Note there are two sets of reads as this is [paired end](https://emea.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) data.


## Preprocessing the reads
### Quality check of the reads
Similarly to what done for the RNA-seq data, the first thing to do in a sequencing project is to check the quality of the data.
We can check the quality of Illumina data is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) it for all our samples
using a `for` loop.
First, we create the output directory named `FASTQC` with `mkdir`:
```
mkdir FASTQC
```

Then, we activate the right anaconda environment:
```
conda activate readmapping_env
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
mkdir REFERENCE
```

To facilitate the process, we provide a reference sequence comprising of only chromosome 20:
```
wget -O REFERENCE/chr20.fa https://www.dropbox.com/s/dq2qgdko6pgk14o/chr20.fa
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
bwa-mem2 mem -R "@RG\tID:Holstein_1\tSM:Holstein_1" -t 1 REFERENCE/chr20.fa TRIM/Holstein_1_R1.trimmed.fq.gz TRIM/Holstein_1_R2.trimmed.fq.gz > ALIGN/Holstein_1.raw.sam
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
fastdup \
   --input FILTER/Holstein_1.filt.bam \
   --output MARKDUP/Holstein_1.md.bam \
   --metrics MARKDUP/Holstein_1.fastdup.stats.txt \
   --create-index \
   --num-threads 2
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
6. [octopus](): a variant caller using tree-based approaches to identify variants.

In our case, we are going to use Clair3, which is a popular, well-established variant caller with plenty of 
functionality. Clair3 has a large number of parameters that can be fine-tuned to achieve greater sensitivity, 
to improve the calling, or to processed pooled samples at once. 

To run Clair3 we need to first download the appropriate model for the illumina reads. The officially trained models can be found [here](https://github.com/HKU-BAL/Clair3#HKU-provided-Models). We need to download the illumina one as follow:
```
mkdir -p models/ilmn
wget -O models/ilmn/pileup.pt https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/ilmn/pileup.pt
wget -O models/ilmn/full_alignment.pt https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/ilmn/full_alignment.pt
```

We then activate the variant calling environment:
```
conda activate variant_calling_env
```

To perform the variant calling we are going to create a new folder where to save our output VCFs:
```
mkdir CLAIR3
```

Then, we can proceed with the variant calling. In our case, we are going to run the caller with base parameters.
We can target the region of interest by providing the region in the BED format:
```
echo -e '20\t31200000\t32800000' > CLAIR3/region.bed
```

We can call the variants as follow: 
```
run_clair3.sh \
  --bam_fn=MARKDUP/Holstein_1.md.bam \
  --ref_fn=REFERENCE/chr20.fa \
  --threads=1 \
  --platform="ilmn" \
  --bed_fn="CLAIR3/region.bed" \
  --model_path="${PWD}/models/ilmn" \
  --sample_name Holstein_1 \
  --gvcf \
  --output=CLAIR3/Holstein_1
```
This will generate an output VCF file ready for downstream analyses.

## Repeat for the second sample
We can now repeat what we have done for a second sample, availble here:
```
wget -O DATA/NDama_1_R1.fq.gz https://www.dropbox.com/s/sl25qu47sghh1zr/NDama_1_R1.fq.gz
wget -O DATA/NDama_1_R2.fq.gz https://www.dropbox.com/s/w1rg72pt12gt2mk/NDama_1_R2.fq.gz
```

Rerun the steps for this second individual, aiming to generate a second GVCF file that we can combine with the previous one.
Remember to activate the right environment to perform the different steps:
```
conda activate readmapping_env
```

and then:
```
conda activate variant_calling_env
```
when performing the variant calling.

> *How many variants does this individual present?*
> *What is the percentage of mapped reads?*
> *If higher/lower than the previous one, what can be causing it?*

## Combining the GVCFs
After variant calling the individual, we can proceed at combining the variants for the two individuals in a single, multi-sample VCF file. We can do so using the [GLNexus](https://github.com/dnanexus-rnd/GLnexus) tool, that provides high speed and compatibility with most tools:
```
mkdir VCF/
```
We then download the appropriate configuration for GLNexus to work with Clair3's outputs:
```
wget -O VCF/clair3.yml http://www.bio8.cs.hku.hk/clair3_trio/config/clair3.yml
```

Before running GLNexus, we need to ensure that each file has a unique name, otherwise the tool will complain:
```
mkdir GVCFS
cp CLAIR3/Holstein_1/merge_output.gvcf.gz GVCFS/Holstein_1.gvcf.gz
cp CLAIR3/Holstein_1/merge_output.gvcf.gz GVCFS/Holstein_1.gvcf.gz.tbi
cp CLAIR3/NDama_1/merge_output.gvcf.gz GVCFS/NDama_1.gvcf.gz
cp CLAIR3/NDama_1/merge_output.gvcf.gz GVCFS/NDama_1.gvcf.gz.tbi
```

and then, we can run glnexus with the desired configuration:
```
glnexus_cli --threads 1 --config VCF/clair3.yml GVCFS/*.g.vcf.gz > VCF/joint_calling.bcf
```
The resulting file will be a single BCF, containg all the variants for all individuals processed and provided to the analysis.
We can visualize its content using [`bcftools view`](https://github.com/samtools/bcftools):
```
bcftools view VCF/joint_calling.bcf | less -S
```

Before starting to work on the file, we will convert it to compressed VCF (`.vcf.gz`) and index it.
It is best practice to compress and index the VCF files generated to save space,
while speeding up the access to the entries of the file. 
We will convert the file to compressed VCF using again [`bcftools view`](https://github.com/samtools/bcftools), while writing the index for it contextually:
```
bcftools view -O z -o VCF/joint_calling.vcf.gz --write-index=tbi VCF/joint_calling.bcf
```

## VCF metrics and statistics
The VCF files generated by any variant caller are large files that store a large amount of information.
This makes them very hard to screen by hand, and not easy to visualise. To have a peek to the file themself, we can use the `less` command:
```
less VCF/joint_calling.vcf.gz
```
Since the lines of the file are long, it can be convenient to add the `-S` option to `less`, which avoid the automatic new line for each row in the VCF:
```
less -S VCF/joint_calling.vcf.gz
```
In this case, we can use all the arrows to move around the file vertically (up/down) and horizontally (right/left).

Given the large size of these files, we can use some software to collect metrics to describe the file.
One of the commonly used tools is [bcftools](https://samtools.github.io/bcftools/bcftools.html), which can screen the VCF for us
and print out generic metrics: 
```
bcftools stats VCF/joint_calling.vcf.gz > VCF/joint_calling.bcftools_stats
```
Remember, the `>` character means "**redirection**": it redirect the output of a command from the screen (known as STDOUT) to a 
file of destination.

> *How many variants does the VCF contains?*
> *How many of each class?*
> *What is the transition/transvertion rate?*


## Filter variants
Combining the files generates a VCF file with a large nubmer fo entries. Several of these might have low coverage issues, low quality or high missingness rate.
These can be filtered out using [vcftools](vcftools.github.io/). We first exclude variants with low QUAL score, suggesting that these are likely wrong calls:
```
vcftools --gzvcf VCF/joint_calling.vcf.gz --minQ 20 --recode --recode-INFO-all --out VCF/joint_calling.highQ
```

Then, we filter out variants with too low depth of sequencing. The depth of sequencing is the number of reads aligned to a given base,
and carrying either the reference or alternative allele. Regions with too low depth of sequencing are more difficult to distinguish
from errors, whereas regions with too high depth are likely to be repetitive regions, enriched with sequencing errors.
The VCF file codify the depth of sequencing in the `DP` field, and is associated to each genotype for each individual.  
To filter out variants with low (depth < 5) or high (depth > 50) depth of sequencing, we use `vcftools` again:
```
vcftools --vcf VCF/joint_calling.highQ.recode.vcf --min-meanDP 5 --max-meanDP 60 --recode --recode-INFO-all --out VCF/joint_calling.highQ.5-60DP
```

Additional possible filtering are accessible in the [vcftools documentation](https://vcftools.github.io/man_latest.html).

> *How many variants do we save after the filtering are applied?*

## Functional variant annotation
Following the quick filtering of the VCF file, we can annotate the variants using either:
1. the Variant Effect Predictor ([VEP](https://www.ensembl.org/info/docs/tools/vep)) tool; or
2. [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction)

`VEP` is a complex software, with multiple options available and configuration that can be used. It also undergoes to frequent revision to incorporate new data, making it important to point out which version of the software we are using. In our case, we are processing the data using the latest version of the software (v115 at the time of writing this tutorial). 
snpEff is also a valid choice, providing easier adoption and good performance while supporting a large collection of pre-compiled annotations.

Then, we annotate the vcf file using the `snpEff` software:
```
snpeff bosTau9 VCF/joint_calling.highQ.5-60DP.recode.vcf -htmlStats VCF/joint_calling.highQ.5-60DP.recode.annot.html -csvStats VCF/joint_calling.highQ.5-60DP.recode.annot.csv > VCF/joint_calling.highQ.5-60DP.recode.annot.vcf
```

This command generate an annotated VCF file.

> *How many variants with moderate effect do we have?*

## Create a summary report with MultiQC
Finally, we can collect generic metrics for the whole analysis using [MultiQC](https://multiqc.info/docs/).
To do so, we can simply run:
```
multiqc .
```
The resulting report will be saved as `multiqc_report.html`, and can be visualized using any browser.