.PHONY : strelka snpeff_download_db 

STRELKA_PATH = $(HOME)/usr/strelka/1.0.15/bin
SNPEFF_PATH = $(HOME)/usr/snpeff/4.3

#----------
# Setup SnpEff
#----------

SNPEFF_GENOME = GRCh37.75

# Download the SnpEff genome database
snpeff_download_db : refs/snpeff/$(SNPEFF_GENOME).timestamp

refs/snpeff/$(SNPEFF_GENOME).timestamp : 
	java -Xmx4G -jar $(SNPEFF_PATH)/snpEff.jar \
		download \
		-c $(SNPEFF_PATH)/snpEff.config \
		$(SNPEFF_GENOME) && touch $@

#----------
# Running Strelka
#----------
strelka : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf 

strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf :
	$(STRELKA_PATH)/configureStrelkaWorkflow.pl \
    --tumor bams/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal bams/HCC1395_exome_normal.17.7MB-8MB.bam \
    --ref ~/share/references/genomes/gsc/GRCh37-lite.fa \
    --config strelka/config/strelka_config_bwa_default.ini \
    --output-dir strelka/HCC1395_exome_tumour_normal

strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.snpeff.vcf : strelka/HCC1395_exome_tumour_normal/results/passed.somatic.snvs.vcf
	java -Xmx4G -jar $(SNPEFF_PATH)/snpEff.jar \
		-canon \
		$(SNPEFF_GENOME) \
		-s $$(basename $$@).summary.html \
		$< > $@

#----------
# Knit Rmarkdown Page
#----------
%.html : %.Rmd
	Rscript --slave -e "rmarkdown::render('$<')"
