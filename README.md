# Variant Calling in Cancer Genomes Seminar

* Date: Sept 27, 2016
* Time: 12-1pm
* Location: Dorothy Lam Boardroom, British Columbia Cancer Research Centre, Vancouver, BC, Canada.

This repository provides instructions on how to perform variant calling in cancer genome data. Specifically, the dataset is a matching tumour and normal exome from a breast cancer cell-line (HCC1395). The data is available from https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data. 

## Contact

Feel free to contact me for help regarding the content in this workshop:

* email: fongchunchan@gmail.com
* twitter: https://twitter.com/fongchunchan
* blog: http://tinyheero.github.io

## Table of Contents

* [Setup](#setup)
    + [Clone Repository](#clone-repository)
    + [Using Conda](#using-conda)
    + [Downloading Human Reference](#downloading-human-reference)
    + [Getting the Full Exome Data](#getting-the-full-exome-data)
        - [Bam to Fastq Conversion](#bam-to-fastq-conversion)
        - [Sequence Alignment using BWA](#sequence-alignment-using-bwa)
        - [Post-Processing the Alignments](#post-processing-the-alignments)
    + [Installing MutationSeq](#installing-mutationseq)
    + [Installing Strelka](#installing-strelka)
    + [Installing SnpEff](#installing-snpeff)
* [Calling Variants](#calling-variants)
    + [Using MutationSeq](#using-mutationseq)
    + [Using Strelka](#using-strelka)
* [Annotating Variants](#annotating-variants)
* [Converting VCF to Table](#converting-vcf-to-table)
* [Post-Processing in R](#post-processing-in-r)
* [Pipeline](#pipeline)

## Setup

These instructions have been tested on a linux machine. 

### Clone Repository

git clone this repository:

```{bash}
git clone git@github.com:tinyheero/variant_calling_in_cancer_genomes_seminar.git
```

The repository provides the following files:

* `bams`: These are the BWA aligned exome bam files that will be used in this tutorial. The bam files have been restricted to a 1 MB region on chromosome 17.
* `Makefile`: A makefile pipeline that executes the commands in this workshop
* `analyze_snv_results.Rmd`: A [rmarkdown](http://rmarkdown.rstudio.com/) file that demonstrates a standard post-processing analysis

### Using Conda

While not necessary, using conda for both installation and package management will make life easier and is recommended for this workshop. Conda is a package management system that is becoming popular in the field of bioinformatics for reproducible research. An increasing number of bioinformatics software are now being distributed through this system. 

For this workshop, installation of different tools will be done by conda (when possible). When the tools are not in conda, instructions on how to manually install the software will be provided. To install conda, we can get it through [miniconda](http://conda.pydata.org/miniconda.html). First download miniconda (for python 2.7) and then run:

```{bash}
sh Miniconda2-latest-Linux-x86_64.sh
```

Then follow the instructions. When you have finished following the instructions, you should have python installed:

```{bash}
which conda
~/miniconda2/bin/conda
```

If you choose not to use conda to install the software needed for this workshop, then you will have to manually install it by yourself. 

### Downloading Human Reference

The first thing we will do is get a human reference to work with. The human reference genome cannot be downloaded from this github repo as the file is too large. But the reference genome can be obtained from the [Genome Science Center FTP server](http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/). We will put these the genome fasta and index file into the refs folder:

```{bash}
mkdir refs
cd refs
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai
```

### Getting the Full Exome Data

> You can skip this section if you are content with working with the bam files that are in the repo. 

The original full exome data can be found https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data. The repo contains in the `bam` folder two smaller tumour and normal bam files where only a 1 MB region on chromosome 17 is represented. This was done to file size issues. If you are interested in working with the whole exome data set, then you can follow these instructions:

```{bash}
cd bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_CGATGT.bam # normal exome
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_ATCACG.bam # tumour exome
```

Once these bam files have been downloaded, you will need to extract them as fastq files. 

#### Bam to Fastq Conversion

Picard SamToFastq provides this functionality. To install this:

```{bash}
conda install -c bioconda picard
picard SamToFastq --version
```

We can now extract to fastq by using the command for the normal exome:

```{bash}
mkdir fastq;
picard SamToFastq \
  INPUT=bam/gerald_C1TD1ACXX_7_CGATGT.bam \
  FASTQ=fastq/gerald_C1TD1ACXX_7_CGATGT_R1.fastq \
  SECOND_END_FASTQ=fastq/gerald_C1TD1ACXX_7_CGATGT_R2.fastq
```

Now for the tumour exome:

```{bash}
picard SamToFastq \
  INPUT=bam/gerald_C1TD1ACXX_7_ATCACG.bam \
  FASTQ=fastq/gerald_C1TD1ACXX_7_ATCACG_R1.fastq \
  SECOND_END_FASTQ=fastq/gerald_C1TD1ACXX_7_ATCACG_R2.fastq
```

These extraction steps will take a fair bit of time. 

#### Sequence Alignment using BWA

Once you have extracted these files, we can align them using bwa. By default bwa outputs sam files and for compression reasons we want to be working in the binary form of sam which is bam. We will also need samtools for this conversion. We can install bwa (v0.7.12) and samtools using conda:

```{bash}
conda install -c bioconda bwa=0.7.12 samtools
```

We will need to first bwa index the genome:

```{bash}
bwa index refs/GRCh37-lite.fa
```

Once you have done this, we can align each exome. For the normal exome:

```
mkdir sai

# bwa aln
bwa aln refs/GRCh37-lite.fa fastq/gerald_C1TD1ACXX_7_CGATGT_R1.fastq > sai/gerald_C1TD1ACXX_7_CGATGT_R1.sai
bwa aln refs/GRCh37-lite.fa fastq/gerald_C1TD1ACXX_7_CGATGT_R2.fastq > sai/gerald_C1TD1ACXX_7_CGATGT_R2.sai

# bwa sampe
bwa sampe refs/GRCh37-lite.fa \
  sai/gerald_C1TD1ACXX_7_CGATGT_R1.sai \
  sai/gerald_C1TD1ACXX_7_CGATGT_R2.sai \
  fastq/gerald_C1TD1ACXX_7_CGATGT_R1.fastq \
  fastq/gerald_C1TD1ACXX_7_CGATGT_R2.fastq | 
  samtools view -bh > HCC1395_exome_normal.bam
```

Now for the the tumour exome:

```
# bwa aln
bwa aln refs/GRCh37-lite.fa fastq/gerald_C1TD1ACXX_7_ATCACG_R1.fastq > sai/gerald_C1TD1ACXX_7_ATCACG_R1.sai
bwa aln refs/GRCh37-lite.fa fastq/gerald_C1TD1ACXX_7_ATCACG_R2.fastq > sai/gerald_C1TD1ACXX_7_ATCACG_R2.sai

# bwa sampe
bwa sampe refs/GRCh37-lite.fa \
  sai/gerald_C1TD1ACXX_7_ATCACG_R1.sai \
  sai/gerald_C1TD1ACXX_7_ATCACG_R2.sai \
  fastq/gerald_C1TD1ACXX_7_ATCACG_R1.fastq \
  fastq/gerald_C1TD1ACXX_7_ATCACG_R2.fastq | 
  samtools view -bh > HCC1395_exome_tumour.bam
```

This would give us the full tumour and normal exome of the HCC1395 cell-line. 

### Post-Processing the Alignments

For this workshop, we will be working with only a 1 MB window of chromosome 17. These bam files are in this repo:

* `bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam`
* `bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam`

These bam files were generated with the following commands. First we sort the samples:

```{bash}
picard SortSam \
  I=bam/HCC1395_exome_normal.bam \
  O=bam/HCC1395_exome_normal.sort.bam \
  SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT 
```

Then we mark for PCR duplicates:

```
picard MarkDuplicates \
  I=bam/HCC1395_exome_normal.sort.bam \
  O=bam/HCC1395_exome_normal.sort.markdup.bam \
  M=bam/markdup_stats/HCC1395_exome_normal_marked_dup_metrics.txt \
  VALIDATION_STRINGENCY=LENIENT 
```

And then finally, we filter for the region of interest:

```{bash}
samtools view -b bam/HCC1395_exome_normal.sort.markdup.bam 17:7000000-8000000 > bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam
```

These commands refer to the normal sample, but are directly applicable to the tumor sample.

### Installing MutationSeq

MutationSeq can be downloaded from http://compbio.bccrc.ca/software/mutationseq. For this workshop, version 4.3.8 was used. Once you have downloaded it, extract it:

```{bash}
tar -xzvf museq_4.3.8.tar.gz
```

This will extract the content into a folder `mutationseq`. We will move the files into `$HOME/usr/museq/4.3.8` for organization purposes:

```
mkdir -p $HOME/usr/museq/4.3.8
mv mutationseq/* $HOME/usr/museq/4.3.8
```

Now we need to install MutationSeq. MutationSeq requires python (v2.7) and several key package dependencies. 

* numpy
* scipy
* matplotlib
* scikit-learn
* intervaltree

The best way to install all of this is to use either [Miniconda](http://conda.pydata.org/miniconda.html) or anaconda. We will use miniconda here. First download miniconda (for python 2.7) and then run:

```{bash}
sh Miniconda2-latest-Linux-x86_64.sh
```

Then follow the instructions. When you have finished following the instructions, you should have python installed:

```{bash}
which python
~/miniconda2/bin/python
```

Now we can install the dependencies needed (NOTE: Versions and specified to ensure compatibility)

```{bash}
conda install -c bioconda numpy=1.7.1 scipy=0.12.0 scikit-learn=0.13.1 matplotlib=1.2.1 intervaltree
```

One last thing that is needed before we can install MutationSeq is the Boost C libraries. We only need to download them from http://www.boost.org/. Once you have downloaded (tested on 1.51) just extract them to a location. For example, you could put it into `$HOME/usr/boost/1.51`

Once this has been installed, we can now proceed to compiling a dependency `pybam.so` (this may be provided in the .tar.gz, but it's better to compile on your own system to ensure that it does work).

```{bash}
make PYTHON=python BOOSTPATH=$HOME/usr/boost/1.51 -B
```

Now when you run:

```{bash}
python $HOME/usr/museq/4.3.8/museq/classify.py --version
4.3.8
```

This indicates that you have successfully installed MutationSeq.

> An important thing to note is that MutationSeq comes bundled with a trained classifer using the [scikit-learn](http://scikit-learn.org/) library and packaged as a pickle. The pickle is only compatible with the specific scikit-learn version that it was built with. As we are using miniconda here with scikit-learn version 0.13.1, we will have to use the corresponding models associated with that version. See [Calling Variants - Using MutationSeq](#using-mutationseq) for more details.

### Installing Strelka

Strelka can be downloaded from https://sites.google.com/site/strelkasomaticvariantcaller/home/download as a .tar.gz file. For this workshop, version 1.0.15 was used. Once you have it downloaded, extract it:

```{bash}
tar -xzvf strelka_workflow-1.0.15.tar.gz
```

This will create a `strelka_workflow-1.0.15` folder. Go into the folder now and run:

```{bash}
cd strelka_workflow-1.0.15
./configure --prefix=$HOME/usr/strelka/1.0.15
make
```

This will install strelka into your home directory folder at `$(HOME)/usr/strelka/1.0.15`. 

Because we are working with exome data in this workshop, we need to make a small adjustment to the Strelka configuration file. After you have finished installing Strelka, there will be a configuration file `strelka_config_bwa_default.ini` (this is Strelka bwa specific configuration). In this file, there is a line:

```
isSkipDepthFilters = 0
```

According to the Strelka FAQ page:

> The depth filter is designed to filter out all variants which are called above a multiple of the mean chromosome depth, the default configuration is set to filter variants with a depth greater than 3x the chromosomal mean. If you are using exome/targeted sequencing data, the depth filter should be turned off…
>
> However in whole exome sequencing data, the mean chromosomal depth will be extremely small, resulting in nearly all variants being (improperly) filtered out.

Sine we are working with exome data here, this value should be set to 1:

```
isSkipDepthFilters = 1
```

The `strelka/config/strelka_config_bwa_default.ini` in this repo is provided to demonstrate this. We will use this file for the workshop.

### Installing SnpEff

SnpEff can be downloaded from https://sourceforge.net/projects/snpeff/files. For this workshop, we will use [version 4.2 of SnpEff](https://sourceforge.net/projects/snpeff/files/snpEff_v4_2_core.zip/download). Once you downloaded unzip the files to

```
unzip snpEff_v4_2_core.zip -d $(HOME)/usr/snpeff/4.3
mv $(HOME)/usr/snpeff/4.3/snpEff/* $(HOME)/usr/snpeff/4.3 # for organization
```

You need to edit the `$(HOME)/usr/snpeff/4.3/snpEff.config` file to add the path where the SnpEff genome databases will be installed. So open this file in your favourite text editor (e.g. vim) and edit the data.dir line to indicate where you want your snpeff databases to be installed.

```
data.dir = /path/to/snpeff/databases
```

For example, we can specify the databases to be installed in a refs folder in your home directory:

```
data.dir = $(HOME)/refs/snpeff/4.3
```

Now we need to build the genome database that we will annotate against. Specifically, we will use the GRCh37.75 genome.

```
java -Xmx4G -jar $(HOME)/usr/snpeff/4.3/snpEff.jar \
  download \
  -c $(HOME)/usr/snpeff/4.3/snpEff.config \
  GRCh37.75
```

## Calling Variants

### Using MutationSeq

```
mkdir -p museq/vcf; \
python $(HOME)/usr/museq/4.3.8/museq/classify.py \
  normal:bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam \
  tumour:bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam \
  reference:refs/GRCh37-lite.fa \
  model:$(HOME)/museq/4.3.8/museq/models_anaconda/model_v4.1.2_anaconda_sk_0.13.1.npz \
  -c $(HOME)/usr/museq/4.3.8/museq/metadata.config \
  -o museq/vcf/HCC1395_exome_tumour_normal_17.vcf
```

### Using Strelka

To call variants in Strelka, we use the following command:

```
$(HOME)/usr/strelka/1.0.15/bin/configureStrelkaWorkflow.pl \
  --tumor bams/HCC1395_exome_tumour.17.7MB-8MB.bam \
  --normal bams/HCC1395_exome_normal.17.7MB-8MB.bam \
  --ref refs/GRCh37-lite.fa \
  --config strelka/config/strelka_config_bwa_default.ini \
  --output-dir strelka/HCC1395_exome_tumour_normal_17
```

This will setup the Strelka run, now do the following:

```
cd strelka/HCC1395_exome_tumour_normal_17
make
```

This will run Strelka and take some time depending on how fast your computer is. 

> You can use the -j parameter to specify how many computational cores to use to speed up the process. For example `make -j 4` will use 4 computational cores. 

## Annotating Variants

We will annotate variants using [SnpEff](http://snpeff.sourceforge.net/). Other options one could try are:

* [Annovar](http://annovar.openbioinformatics.org/en/latest/)
* [VEP](http://uswest.ensembl.org/info/docs/tools/vep/index.html)

Once we have the vcf files from MutationSeq and Strelka, we can use the following commands:

```
# Annotating MutationSeq Vcf
java -Xmx4G -jar $(HOME)/usr/snpeff/4.3/snpEff.jar \
  -canon \
  GRCh37.75 \
  -s museq/vcf/HCC1395_exome_tumour_normal.snpeff.summary.html \
  museq/vcf/HCC1395_exome_tumour_normal_17.vcf \
  > museq/vcf/HCC1395_exome_tumour_normal.snpeff.vcf

# Annotating Strelka Vcf
java -Xmx4G -jar $(HOME)/usr/snpeff/4.3/snpEff.jar \
  -canon \
  GRCh37.75 \
  -s strelka/HCC1395_exome_tumour_normal_17/results/passed.somatic.snvs.snpeff.summary.html \
  strelka/HCC1395_exome_tumour_normal_17/results/passed.somatic.snvs.vcf \
  > strelka/HCC1395_exome_tumour_normal_17/results/passed.somatic.snvs.snpeff.vcf
```

## Converting VCF to Table

One standard step that is often done is converting the VCF file into a tabular format that can be easily loaded into other software (e.g. R, python) for additional analysis or easier distribution (e.g. excel) for collaborators. 

The following command demonstrates how one can convert the Strelka VCF output file that has been annotated with SnpEff into a tabular format.

```
java -jar $(HOME)/usr/snpeff/4.3/SnpSift.jar \
  extractFields \
  -e "."  \
  -s "," \
  strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.vcf \
  CHROM POS ID REF ALT QUAL FILTER QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC \
  GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP \
  GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU \
  GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU ANN[*].ALLELE ANN[*].EFFECT \
  ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID \
  ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS \
  ANN[*].CDNA_LEN ANN[*].CDS_POS ANN[*].CDS_LEN ANN[*].AA_POS \
  ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS \
  > strelka/HCC1395_exome_tumour_normal_17/results/passed.somatic.snvs.snpeff.tsv
```

## Post-Processing in R

The final step is often the post-processing of the results in a data analysis language. In this workshop, we will use the data analysis language R for our post-processing. The files we will be post-processing are in this repo:

* `strelka/HCC1395.strelka.full.txt`
* `museq/HCC1395.museq.full.txt`

These are the Strelka and MutationSeq runs on the full exome data as opposed to the subset of the exome which are what the bam files in the repository are. The full exome data though is processed through the same pipeline though.

A rmarkdown file (`analyze_snv_results.Rmd`) has been provided that provides the R code to demonstrate some of the typical plots and analyses that can be generated from variant calling results. The rmarkdown file can be rendered into a html page that can be opened in a standard web browser (e.g. Google Chrome). You will need Rstudio (v0.99.903; tested on this version) in order to render the rmarkdown file. Also the following R packages need to be installed:

* data.table (v1.9.6)
* ggplot2 (v2.1.0)
* plyr (v1.8.3)
* dplyr (v0.5.0)
* stringr (v0.6.2)

Now when you open the `analyze_snv_results.Rmd` file in Rstudio, you should be able to press the "Knit HTML") button and it should render the rmarkdown file into a html page.

## Pipeline

All of these steps have been packaged into a Makefile pipeline that is in this github repo. You just need to complete the "Setup" section and then you should be able to just run:

```
make variant_call
```

And this will perform the variant calling and annotation steps for you.
