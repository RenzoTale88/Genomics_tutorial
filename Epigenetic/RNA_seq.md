# RNA-seq analysis tutorial
Adapted from James Prendergast's lecture at the 2022 GCRF-STAR workshop in Arusha, Tanzania.

In this tutorial we will go over how to analyse Illumina RNA-seq data. This will comprise four main steps:
1. Sequencing data QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
1. Alignment of sequencing data to genome ([STAR](https://github.com/alexdobin/STAR))
1. Gene level quantification ([SALMON](https://github.com/COMBINE-lab/salmon))
1. Downstream analysis such as gene differential expression ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))

## Setting up the environment

To install the various tools we will need to analyse the RNA-seq data we will need to first install a range of packages using Ubuntu’s package manager apt. These are required for, for example, to compile the software. On other servers many of these packages may already be installed. Once more, we will do it using [anaconda](https://anaconda.org/).

To install mamba, please follow the instructions [here](../VariantCalling/WGS_tutorial.md).

Then, run the virtual environment creation:
```
mamba create -n rnaseq_env -y -c conda-forge -c bioconda samtools multiqc STAR salmon bioconductor-deseq2 bioconductor-pcaexplorer gffread bioconductor-reportingtools bioconductor-tximport bioconductor-genomeinfodbdata r-markdown fastqc
```
And then activate the environment:
```
mamba activate rnaseq_env
```

## Obtaining the sequencing data

For this tutorial we have made available sequencing reads from B cells of a Holstein-Friesian cow on dropbox. We will first make a folder for us to work in (mkdir) then move into this folder (cd):
```
mkdir RNASeq && cd RNASeq
```

We can now download the sequencing data into the new folder. To do this we will use the wget program that can pull files from the internet (in this case dropbox):
```
mkdir -p DATA
wget -O DATA/HF3471_bcell_chr1_1.fq.gz https://www.dropbox.com/s/scwahcxvw92iwhw/HF3471_bcell_chr1_1.fq.gz
wget -O DATA/HF3471_bcell_chr1_2.fq.gz https://www.dropbox.com/s/0o93kgjh0hidkw1/HF3471_bcell_chr1_2.fq.gz
```

Note there are two sets of reads as this is paired end data.

You can view the fastq files using `less`, which is a special version of the `less` tool for viewing compressed files:
```
less HF3471_bcell_chr1_1.fq.gz
```

As with `less` just need to press q to quit the viewer.

## Data quality control

The first thing to do for almost any sequencing project is to check the quality of the data, as, for example, you may need to do some filtering/trimming if find any problems with it. A widely used tool for checking the quality of Illumina data is FastQC.

To run FastQC on the sequencing data we downloaded from Dropbox we first make a directory to store the results in (using mkdir) then run FastQC specifying the sequencing files and this output directory:
```
mkdir FastQC_output
fastqc DATA/HF3471_bcell_chr1_1.fq.gz DATA/HF3471_bcell_chr1_2.fq.gz --outdir=FastQC_output/
```

There should now be files ending in html in the FastQC_output folder.

Have a look at your FastQC reports.

> *How many sequence reads are there in each file?*

> *Which metrics do you have warnings on?*

## Prepare for alignment

To align our paired-end reads we will use the (STAR) aligner that is designed specifically for aligning RNA-seq data. 
In our case, we will use the STAR aligner. 

To run STAR we first need to make a genome index. To do this we need to download the genome we want to align against and a gtf file specifying the location of known genes. STAR uses the latter to determine where known splice sites are (though STAR is still able to find other splice sites in the reads).

To make things run more quickly, rather than align reads to the whole cow genome, we will just map them to chromosome 1. We have already extracted the sequence of chromosome 1 and saved it in the chr1.fa file on dropbox.
```
mkdir REFERENCE
wget -O REFERENCE/chr1.fa.gz https://www.dropbox.com/s/sn5grtpczlgf48l/chr1.fa.gz
```
We will download the gtf file specifying the location of known cow genes from Ensembl and then uncompress both files.
```
wget -O REFERENCE/Bos_taurus.ARS-UCD1.2.106.gtf.gz https://www.dropbox.com/s/1i145zewlw43gjn/Bos_taurus.ARS-UCD1.2.106.gtf.gz
gunzip REFERENCE/chr1.fa.gz
gunzip REFERENCE/Bos_taurus.ARS-UCD1.2.106.gtf.gz
```
Take a look at these files to see their format (remember to press q to exit the viewer):
```
less REFERENCE/chr1.fa
less REFERENCE/Bos_taurus.ARS-UCD1.2.106.gtf
```

More details on the format of the gtf file can be found here.

We will now use these to make a genome index to align our reads against. First we will make a directory to store the index in:
```
mkdir STAR_genome
```

And then we build the index:
```
STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir STAR_genome \
--genomeSAindexNbases 12 \
--genomeFastaFiles REFERENCE/chr1.fa \
--sjdbGTFfile REFERENCE/Bos_taurus.ARS-UCD1.2.106.gtf \
--sjdbOverhang 99
```

Can see we specified a range of different parameters alongside specifying the two files we just downloaded. An important one is `–runMode genomeGenerate` which tells STAR we want to build an index. Have a look at the manual to see what each of these parameters are for.

> *What does specifying `–runThreadN 3` do?*
> *Align the reads*

Now we have the genome index we can align our RNA-seq reads to it. As before we first make a folder to store the output in with mkdir, then run the alignment with STAR:
```
mkdir STAR_output

STAR --runThreadN 3 \
--genomeDir STAR_genome \
--readFilesIn DATA/HF3471_bcell_chr1_1.fq.gz DATA/HF3471_bcell_chr1_2.fq.gz \
--readFilesCommand zcat \
--quantMode TranscriptomeSAM GeneCounts \
--outFileNamePrefix ./STAR_output/
```

You can check the STAR_output folder for the STAR results.
```
ls STAR_output
```

There are a number of files in here. Sam and bam files containing the read alignments, read counts by genes and log files summarising the results. Lets have a look at one of the log files using the `less` command:
```
less STAR_output/Log.final.out
```
To scroll down the file press space and to quit `less` press `q`.

> *How many reads were successfully mapped?*
> *Quanitfy gene expression*

Now we have successfully aligned our genes to the genome and transcriptome we are going to quantify each gene’s expression level. There are many ways to quanitfy gene expression (including within STAR) but for this tutorial we will use Salmon.

Salmon needs the sequences of the transcripts used. We will extract these from the gtf and genome fasta file we used with STAR. Using the `gffread` software we extract the sequences for the transcripts. In this tutorial we are only using genes on chromosome 1 so first we pull out of the gtf file those genes on this chromosome. Then use gffread to get the sequences:

```
grep -w "^1" REFERENCE/Bos_taurus.ARS-UCD1.2.106.gtf > REFERENCE/chr1.gtf
gffread -w REFERENCE/transcripts.fa -g REFERENCE/chr1.fa REFERENCE/chr1.gtf
```

Now if you look at the `transcripts.fa` file can see contains the sequences of the transcripts:
```
less REFERENCE/transcripts.fa
```
Now we can use this with the output from STAR to get transcript level expression levels
```
salmon quant --threads 3 \
        --targets REFERENCE/transcripts.fa \
        --libType A \
        --output salmon_output/HF3471_bcell \
        --alignments ./STAR_output/Aligned.toTranscriptome.out.bam
```

Salmon should now have quantified the expression levels of each transcript. You can see the expression level estimates in the file with the `.sf` extension:
```
less salmon_output/HF3471_bcell/quant.sf
```

## Summarising results with MultiQC

[MultiQC](https://multiqc.info/) is a software tool that can summarise across results from various tools we used in this tutorial. We can just run it in our folder:
```
multiqc .
```
Amd the tool will search through the directory and sub-directories to find the output of the different tools and generate a report, `multiqc_report.html`, that you can transfer to your computer and open.

Now you have successfully processed the RNA-seq reads by checking their quality, mapping them to transcripts and quantifying them. Normally you would have more than one sample and repeat this process for each. In the next section we will show you how to compare the salmon output of groups of samples to identify which genes are differentially expressed between them.
Differential expression tutorial

## Getting data and setting up R

To test for differential expression of genes between groups of samples we are going to use the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package in [R](https://www.r-project.org/). As in the previous section we only processed one sample, we are going to delete the old salmon results folder and replace it with results for six samples that were processed using the same pipeline but with the whole genome rather than just chromosome 1.
```
rm -R salmon_output/
wget https://www.dropbox.com/s/kyah6svaygsnuuc/salmon_output.tar.gz
tar -xzvf salmon_output.tar.gz
```

If we look in this directory you can see we now have six directories containing salmon output for three Nelore cattle and three Holstein:
```
ls salmon_output
ls salmon_output/Nelore1/
```

Then, we can start R by typing:
```
R
```
Within R we first need to install the R packages we are going to use from Bioconductor, including DESeq2. The RNA-seq data analysis environment comes with all the dependencies to run a differential expression analysis with DESeq2. We can simply load the packages using the `library` command:
```
library(biomaRt)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(ReportingTools)
```
We are now ready to go!

### Getting differentially expressed genes with DESeq2

When we ran salmon it quantified the expression level of each transcript. However, DESeq2 requires us to provide information on which gene each transcript belongs to. To get this information we will query the Ensembl database using the biomaRt package. It requires us to specify which dataset at Ensembl we want to query. As this is cows we specify btaurus_gene_ensembl, but for other species you would change accordingly.
```
txdb1 <- txdbmaker::makeTxDbFromBiomart(dataset="btaurus_gene_ensembl")
txdb1
```

The txdbl object now contains information on the ensembl transcripts we used. We can extract the transcript to gene mappings from this using the following code:
```
k <- keys(txdb1, keytype = "TXNAME")
tx2gene <- select(txdb1, k, "GENEID", "TXNAME")
```
We can view the top of this tx2gene object using the following code:
```
head(tx2gene)
```
The column on the left is the Ensembl transcript ids and the column on the right the corresponding gene ids.

We now need to find the location of the salmon output files containing the read counts for each transcript. Rather than specify the locations of these files manually we can use the list.files function to search for them within the current directory:
```
files<-list.files(pattern = ".sf$", recursive = TRUE)
```
If we then view the content of the files object it should show the location of the six salmon output files
```
files
```
We can then label each of the files with the name of the animal and read the files in:
```
names(files)<-c("Holstein1", "Holstein2", "Holstein3", "Nelore1", "Nelore2", "Nelore3")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
head(txi$abundance)
```
You should have seen the head command has printed the expression levels of the first six genes in each animal.

We are now going to generate a table specifying which animal is from which breed:
```
samps<- names(files)
breed <- factor(rep(c("Holstein","Nelore"),each=3))
samples<-data.frame(samps, breed)
rownames(samples) <- samples$samps
samples
```
Now we can use this information to run the differential expression analysis between the two breeds, to identify which genes show the biggest expression difference between Holstein and Nelore. First we read the data into a DESeq2 object.
```
dds <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ breed)
```
The “design = ~ breed” line means we want to make a comparison between the groups specified in the breed column of the samples table we made above.

We will first exclude lowly expressed genes (those with less than ten reads across samples) to increase speed and as we will not be able to reliably detect differential expression at these genes:
```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
We then run the differential expression analysis:
```
dds <- DESeq(dds)
res <- results(dds)
res
```
To make the results easier to view we can finally generate a report listing the most differentially expressed genes:
```
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory = "./DESeq2_output")
publish(dds,des2Report, pvalueCutoff=0.00001, factor = colData(dds)$breed, reportDir="./DESeq2_output")
finish(des2Report)
```
A directory should have been generated, DESeq2_output, that contains an html file reporting the most significant genes. To view this copy the folder over to your computer using cyberduck and open the html report. The file should look something like this:

Thats it! So have gone from raw sequencing reads to a list of differentially expressed genes.

### Sample QC

If have time you can try plotting a heatmap and PCA of the relationships between the samples. You can do this with the PCAexplorer package however to run this you will need to do it on your laptop rather than the server, as it generates an interactive panel. So would download the example salmon data from dropbox to your laptop and rerun the R code as above. Then can run the code below to generate the report:
```
library("pcaExplorer")
pcaExplorer(dds = dds) %>% ggsave("rnaseq_pca.pdf")
```
