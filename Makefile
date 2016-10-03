.PHONY : bwa.index snpeff_download_db fastq bam bam.markdup strelka 

STRELKA_PATH = $(HOME)/usr/strelka/1.0.15/bin
SNPEFF_PATH = $(HOME)/usr/snpeff/4.3

MUSEQ_PATH = $(HOME)/usr/museq/4.3.8/museq

#----------
# Setup Reference Genome
#----------
bwa.index :
	bwa index refs/GRCh37-lite.fa 

#----------
# Setup SnpEff
#----------

SNPEFF_GENOME = GRCh37.75

# Download the SnpEff genome database
snpeff_download_db : snpeff_$(SNPEFF_GENOME).timestamp

snpeff_$(SNPEFF_GENOME).timestamp : 
	java -Xmx4G -jar $(SNPEFF_PATH)/snpEff.jar \
		download \
		-c $(SNPEFF_PATH)/snpEff.config \
		$(SNPEFF_GENOME) && touch $@

#----------
# Download Full Exome Data
#----------
download : bam/gerald_C1TD1ACXX_7_CGATGT.bam bam/gerald_C1TD1ACXX_7_ATCACG.bam

# Exome Normal
bam/gerald_C1TD1ACXX_7_CGATGT.bam :
	mkdir -p $(@D); \
	wget -p $(@D) https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_CGATGT.bam

# Exome Tumor
bam/gerald_C1TD1ACXX_7_ATCACG.bam :
	mkdir -p $(@D); \
	wget -p $(@D) https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_ATCACG.bam

#---------
# Fastq Extraction
#---------
fastq : fastq/gerald_C1TD1ACXX_7_CGATGT_R1.fastq \
	fastq/gerald_C1TD1ACXX_7_CGATGT_R2.fastq \
	fastq/gerald_C1TD1ACXX_7_ATCACG_R1.fastq \
	fastq/gerald_C1TD1ACXX_7_ATCACG_R2.fastq

fastq/%_R1.fastq fastq/%_R2.fastq : bam/%.bam
	picard SamToFastq \
    INPUT=$< \
    FASTQ=fastq/$*_R1.fastq \
    SECOND_END_FASTQ=fastq/$*_R2.fastq

#----------
# Sequence Alignment using BWA
#----------
sai/%.sai: fastq/%.fastq
	mkdir -p $(@D); \
	bwa aln -t 4 refs/GRCh37-lite.fa $< > $@ 2> log/$(@F).log

bam : bam/HCC1395_exome_normal.bam bam/HCC1395_exome_tumour.bam

bam/HCC1395_exome_normal.bam : sai/gerald_C1TD1ACXX_7_CGATGT_R1.sai sai/gerald_C1TD1ACXX_7_CGATGT_R2.sai fastq/gerald_C1TD1ACXX_7_CGATGT_R1.fastq fastq/gerald_C1TD1ACXX_7_CGATGT_R2.fastq
	bwa sampe refs/GRCh37-lite.fa $^ | samtools view -bh > $@

bam/HCC1395_exome_tumour.bam : sai/gerald_C1TD1ACXX_7_ATCACG_R1.sai sai/gerald_C1TD1ACXX_7_ATCACG_R2.sai fastq/gerald_C1TD1ACXX_7_ATCACG_R1.fastq fastq/gerald_C1TD1ACXX_7_ATCACG_R2.fastq
	bwa sampe refs/GRCh37-lite.fa $^ | samtools view -bh > $@

#----------
# Post-Processing of BAM Files
#----------
bam.sort.markdup : bam/HCC1395_exome_normal.sort.markdup.bam bam/HCC1395_exome_tumour.sort.markdup.bam

bam/%.sort.bam : bam/%.bam
	picard SortSam \
		I=$< \
		O=$@ \
		SORT_ORDER=coordinate \
		VALIDATION_STRINGENCY=LENIENT 

bam/%.markdup.bam : bam/%.bam
	mkdir -p bam/markdup_stats; \
	picard MarkDuplicates \
		I=$< \
		O=$@ \
		M=bam/markdup_stats/$*_marked_dup_metrics.txt \
		VALIDATION_STRINGENCY=LENIENT 

# Generate Samtools Index
%.bam.bai : %.bam
	samtools index $<

#----------
# Filter Bam for Only 1MB region of Chromosome 17 
#----------
filtered.bam : \
	bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam \
	bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam.bai \
	bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam \
	bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam.bai

bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam : bam/HCC1395_exome_normal.sort.markdup.bam
	samtools view -b $< 17:7000000-8000000 > $@

bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam : bam/HCC1395_exome_tumour.sort.markdup.bam
	samtools view -b $< 17:7000000-8000000 > $@

#----------
# Variant Calling 
#----------

# Run MutationSeq on Full Exome
museq/vcf/HCC1395_exome_tumour_normal.vcf : bam/HCC1395_exome_normal.sort.markdup.bam bam/HCC1395_exome_tumour.sort.markdup.bam bam/HCC1395_exome_normal.sort.markdup.bam.bai bam/HCC1395_exome_tumour.sort.markdup.bam.bai
	mkdir -p $(@D); \
	python $(MUSEQ_PATH)/classify.py \
		normal:$< \
		tumour:$(word 2,$^) \
		reference:refs/GRCh37-lite.fa \
		model:$(MUSEQ_PATH)/models_anaconda/model_v4.1.2_anaconda_sk_0.13.1.npz \
		-c $(MUSEQ_PATH)/metadata.config \
		-o $@

# Run MutationSeq on Sampled Exome
museq/vcf/HCC1395_exome_tumour_normal_17.vcf : bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam
	mkdir -p $(@D); \
	python $(MUSEQ_PATH)/classify.py \
		normal:$< \
		tumour:$(word 2,$^) \
		reference:refs/GRCh37-lite.fa \
		model:$(MUSEQ_PATH)/models_anaconda/model_v4.1.2_anaconda_sk_0.13.1.npz \
		-c $(MUSEQ_PATH)/metadata.config \
		-o $@

strelka/HCC1395_exome_tumour_normal_17/results/passed.somatic.snvs.vcf : bam/HCC1395_exome_normal.sort.markdup.17.7MB-8MB.bam bam/HCC1395_exome_tumour.sort.markdup.17.7MB-8MB.bam
	$(STRELKA_PATH)/configureStrelkaWorkflow.pl \
    --normal $< \
    --tumor $(word 2,$^) \
    --ref refs/GRCh37-lite.fa \
    --config strelka/config/strelka_config_bwa_default.ini \
    --output-dir strelka/HCC1395_exome_tumour_normal_17

strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf : bam/HCC1395_exome_normal.sort.markdup.bam bam/HCC1395_exome_tumour.sort.markdup.bam
	$(STRELKA_PATH)/configureStrelkaWorkflow.pl \
    --normal $< \
    --tumor $(word 2,$^) \
    --ref refs/GRCh37-lite.fa \
    --config strelka/config/strelka_config_bwa_default.ini \
    --output-dir strelka/HCC1395_exome_tumour_normal && \
	make -C strelka/HCC1395_exome_tumour_normal

#----------
# Annotating Variants
#----------

# Annotate Strelka Results
strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.vcf : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf
	java -Xmx4G -jar $(SNPEFF_PATH)/snpEff.jar \
		-canon \
		-no-downstream \
		-no-intergenic \
		-no-upstream \
		$(SNPEFF_GENOME) \
		-s $(basename $@).summary.html \
		$< > $@

# Convert Strelka VCF to Table
STRELKA_VCF_COLS = CHROM \
									 POS \
									 ID \
									 REF \
									 ALT \
									 QUAL \
									 FILTER \
									 QSS \
									 TQSS \
									 NT \
									 QSS_NT \
									 TQSS_NT \
									 SGT \
									 SOMATIC

SNPEFF_EFF_COLS = GEN[0].DP \
									GEN[1].DP \
									GEN[0].FDP \
									GEN[1].FDP \
									GEN[0].SDP \
									GEN[1].SDP \
									GEN[0].SUBDP \
									GEN[1].SUBDP \
									GEN[0].AU \
									GEN[1].AU \
									GEN[0].CU \
									GEN[1].CU \
									GEN[0].GU \
									GEN[1].GU \
									GEN[0].TU \
									GEN[1].TU
									 
SNPEFF_ANN_COLS = ANN[*].ALLELE \
									ANN[*].EFFECT \
									ANN[*].IMPACT \
									ANN[*].GENE \
									ANN[*].GENEID \
									ANN[*].FEATURE \
									ANN[*].FEATUREID \
									ANN[*].BIOTYPE \
									ANN[*].RANK \
									ANN[*].HGVS_C \
									ANN[*].HGVS_P \
									ANN[*].CDNA_POS \
									ANN[*].CDNA_LEN \
									ANN[*].CDS_POS \
									ANN[*].CDS_LEN \
									ANN[*].AA_POS \
									ANN[*].AA_LEN \
									ANN[*].DISTANCE \
									ANN[*].ERRORS

# -e is the empty field value
# -s multiple value separator
# See http://snpeff.sourceforge.net/SnpSift.html#Extract for more details.
strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.tsv : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.vcf 
	java -jar $(SNPEFF_PATH)/SnpSift.jar \
		extractFields \
		-e "."  \
		-s "," \
		$< $(STRELKA_VCF_COLS) $(SNPEFF_EFF_COLS) $(SNPEFF_ANN_COLS) \
		> $@

#----------
# Run Pipeline
#----------
variant_call : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.tsv

#----------
# Knit Rmarkdown Page
#----------
%.html : %.Rmd
	Rscript --slave -e "rmarkdown::render('$<')"
