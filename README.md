# Variant Calling in Cancer Genomes Seminar

* Date: Sept 27, 2016
* Time: 12-1pm
* Location: Dorothy Lam Boardroom

## Setup

These instructions have been tested on a linux machine. 

### Installing Strelka

Strelka can be downloaded from https://sites.google.com/site/strelkasomaticvariantcaller/home/download as a .tar.gz file. For this workshop, version 1.0.15 was used. Once you have it downloaded, extract it:

```
tar -xzvf strelka_workflow-1.0.15.tar.gz
```

This will create a `strelka_workflow-1.0.15` folder. Go into the folder now and run:

```
cd strelka_workflow-1.0.15
./configure --prefix=$HOME/usr/strelka/1.0.15
make
```

This will install strelka into your home directory folder at `$(HOME)/usr/strelka/1.0.15`. 

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

### Downloading Human Reference

The human reference genome cannot be downloaded from this github repo as the file is too large. But the reference genome can be obtained from the [Genome Science Center FTP server](http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/). We will put these the genome fasta and index file into the refs folder:

```
mkdir refs
cd refs
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai
```

## Calling Variants

### Using Strelka

To call variants in Strelka, we use the following command:

```
STRELKA_PATH=${HOME}/usr/strelka/1.0.15/bin
${STRELKA_PATH}/configureStrelkaWorkflow.pl \
  --tumor bams/HCC1395_exome_tumour.17.7MB-8MB.bam \
  --normal bams/HCC1395_exome_normal.17.7MB-8MB.bam \
  --ref refs/GRCh37-lite.fa \
  --config strelka/config/strelka_config_bwa_default.ini \
  --output-dir strelka/HCC1395_exome_tumour_normal
```

This will setup the Strelka run, now do the following:

```
cd strelka/HCC1395_exome_tumour_normal
make
```

This will run Strelka and take some time depending on how fast your computer is. 

> You can use the -j parameter to specify how many computational cores to use to speed up the process. For example `make -j 4` will use 4 computational cores. 

## Annotating Variants

We will annotate variants using [SnpEff](http://snpeff.sourceforge.net/). Other options one could try are:

* [Annovar](http://annovar.openbioinformatics.org/en/latest/)
* [VEP](http://uswest.ensembl.org/info/docs/tools/vep/index.html)

Once we have the 

## Pipeline

All of these steps have been packaged into a Makefile pipeline that is in this github repo. You just need to complete the "Setup" section and then you should be able to just run:

```
make variant_call
```

And this will perform the variant calling and annotation steps for you.
