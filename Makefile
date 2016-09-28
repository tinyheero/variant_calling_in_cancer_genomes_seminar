.PHONY : strelka snpeff_download_db 

STRELKA_PATH = $(HOME)/usr/strelka/1.0.15/bin
SNPEFF_PATH = $(HOME)/usr/snpeff/4.3

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
# Variant Calling
#----------

# Run Strelka
strelka : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf 

strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf :
	$(STRELKA_PATH)/configureStrelkaWorkflow.pl \
    --tumor bams/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal bams/HCC1395_exome_normal.17.7MB-8MB.bam \
    --ref ~/share/references/genomes/gsc/GRCh37-lite.fa \
    --config strelka/config/strelka_config_bwa_default.ini \
    --output-dir strelka/HCC1395_exome_tumour_normal

# Run MutationSeq

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
# Knit Rmarkdown Page
#----------
%.html : %.Rmd
	Rscript --slave -e "rmarkdown::render('$<')"
