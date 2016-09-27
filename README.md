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

To call variants in strelka, we use the following command:

