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

This will install strelka into your home directory folder at `$HOME/usr/strelka/1.0.15`. 

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
STRELKA_PATH = ${HOME}/usr/strelka/1.0.15/bin
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

