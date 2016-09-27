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
