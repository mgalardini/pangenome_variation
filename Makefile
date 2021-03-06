# Input files
# The genome file should be put in the same directory as the makefile
GENOME = genome.fasta
GBK = genome.gbk
# Scripts directory
SRCDIR = $(CURDIR)/src
# Fles directory
FILESDIR = $(CURDIR)/files
# Directory containing the target genomes fasta files
TARGETSDIR = $(CURDIR)/genomes
# gkm-SVM directory
GKDIR = $(CURDIR)/gkmsvm
# GFF3 directory
GFFDIR = $(CURDIR)/gff
# Output directory for kSNP
KOUT = $(CURDIR)/kout
# Output directory for parsnp
POUT = $(CURDIR)/pout
# Output directory for parsnp (all vs all)
TOUT = $(CURDIR)/tout
# Output directory for reads alignment
MOUT = $(CURDIR)/mout
# Reports dir
NOTEBOOKDIR = $(CURDIR)/notebooks
# Directory where the prokka is
PROKKA = $(SOFTDIR)/prokka-1.10/bin
# Directory where the parsnp binary is
PARSNP = $(SOFTDIR)/Parsnp-Linux64-v1.2
# Directory where the jellyfish binary is
JELLYFISHDIR = $(SOFTDIR)/Jellyfish/bin
# Directory where PICARD is
PICARDDIR = $(SOFTDIR)/picard-tools-1.119
# Directory for GATK
GATKDIR = $(SOFTDIR)/GATK
# Location of JAVA version 7
JAVA7 = $(SOFTDIR)/jre1.7.0_76/bin/java
# How much RAM for JAVA (in Gb)
JAVAMEM = 24
# Directory containing the genomes with paired-ends reads
# Organized in sub-directories
READSDIR = $(CURDIR)/reads
READ1 = READ1.txt.gz
READ2 = READ2.txt.gz
# Proteome fasta files directory
PROTEOMEDIR = $(CURDIR)/proteomes
# Output dir for LS-BSR on reference genes
REFERENCEDIR = $(CURDIR)/lreference
# Output dir for LS-BSR on all genes
ALLDIR = $(CURDIR)/lall
# Directory where the usearch binary is
USEARCHDIR = $(SOFTDIR)/usearch
# Output dir for pairwise oma runs
OMAOUT = $(CURDIR)/oout
# Output dir for kmer counting
KMOUT = $(CURDIR)/kmout
KMLITEOUT = $(CURDIR)/kmliteout
# Output dir for roary
RDIR = $(CURDIR)/roary
# Output dir for roary
SNPEFFDIR = $(CURDIR)/snpEff
# Outgroups (genomes that belong to another species or are very diverse)
OUTGROUPS = outgroups.txt
# Evolution experiments
EVOLUTION = evolution_experiment.txt
# Mash
MASHDIR = $(CURDIR)/mash

# Useful handle for submission to a cluster
SUBMIT = eval

# Parameters
# snpEff annotation name
SNPEFFCHROM = NC_000913.3
# Species (from prokka reannotation)
GENUS = Escherichia
SPECIES = coli
# kSNP CPUs
KCPU = 20
# NOTE: optimal k-mer shoudld be derived from
# a Kchooser run
KKMER = 19
# K-mer counting parameters
WORD = 40
KHASH = 6
# parsnp CPUs
PCPU = 1
# Reads mapping CPUs
MCPU = 2
# LS-BSR CPUS
LCPU = 20
# ROARY CPUS
RCPU = 30
# OMA CPUS
OCPU = 10
# Prokka CPUS
PROKKACPU = 5
# Maximum coverage (for downsampling)
MAXCOVERAGE = 100
SEED = 100
# VCF filters (after reads mapping)
FILTER = -f "DP > 4" -g "GQ > 20"
# Species
SPECIES = ecoli
# Theta (estimated mutation rate)
THETA = 0.05
# How many chromosomal copies?
PLOIDY = 1
# MASH CPUS
MCPU = 5

# Anything below this point should not be changed
$(TARGETSDIR): $(GENOME)
	mkdir -p $(TARGETSDIR)

##################################
## k-mer variant calling (kSNP) ##
##################################

KINPUT = kinput
$(KINPUT): $(GENOME) $(TARGETSDIR)
	echo -e $(CURDIR)/$(GENOME)"\t"$(shell basename $(GENOME)) > $(KINPUT)
	for genome in $$(ls $(TARGETSDIR)); do \
		echo -e $(TARGETSDIR)/$$genome"\t"$$(basename $$genome) >> $(KINPUT); \
	done

KANNOTATED = kannotated
$(KANNOTATED):
	echo $(shell basename $(GENOME)) > $(KANNOTATED)

$(KOUT):
	mkdir -p $(KOUT)

KOUTPUT = $(KOUT)/SNPs_all_matrix
$(KOUTPUT): $(GENOME) $(GBK) $(TARGETSDIR) $(KINPUT) $(KANNOTATED) $(KOUT)
	kSNP3 -vcf -in $(KINPUT) -outdir $(KOUT) -k $(KKMER) -CPU $(KCPU) -annotate $(KANNOTATED) -genbank $(GBK)

KVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(KOUT)/,$(addsuffix .vcf,$(notdir $(basename $(GENOME))))))
KNONSYNVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(KOUT)/,$(addsuffix .nonsyn.vcf,$(notdir $(basename $(GENOME))))))
KTFBSVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(KOUT)/,$(addsuffix .tfbs.vcf,$(notdir $(basename $(GENOME))))))

$(KOUT)/%.vcf: $(TARGETSDIR)/%.fasta $(KOUTPUT)
	$(SRCDIR)/ksnp2vcf $(KOUT)/VCF.$(notdir $(GENOME)).vcf $(notdir $<) --chrom $(shell grep VERSION $(GBK) | head -n 1 | awk '{print $$2}') > $@

$(KOUT)/%.nonsyn.vcf: $(KOUT)/%.vcf $(GBK)
	cat $< | python2 $(SRCDIR)/vcf2nonsyn $(GBK) - > $@

$(KOUT)/%.tfbs.vcf: $(KOUT)/%.vcf $(TFBSTABLE)
	cat $< | python2 $(SRCDIR)/vcf2tfbs $(TFBSTABLE) $(FILESDIR)/pssm - > $@

#################################################
## Alignment variant calling (pairwise parsnp) ##
#################################################

GENOMES = $(wildcard $(TARGETSDIR)/*.fasta)

$(POUT):
	mkdir -p $(POUT)

# Mask repeats in the reference genome
REPEATS = repeats.bed
$(REPEATS): $(GENOME)
	nucmer --maxmatch --nosimplify $(GENOME) $(GENOME) && \
	show-coords -r -T out.delta -H | tail -n+2 > repeats.txt && \
	awk '{print $$8"\t"$$1"\t"$$2}' repeats.txt > $(REPEATS)

PVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .vcf,$(notdir $(basename $(GENOME))))))
PMERGEDVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .merged.vcf,$(notdir $(basename $(GENOME))))))
PANNOTATEDVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .annotated.vcf,$(notdir $(basename $(GENOME))))))
PNONSYNVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .nonsyn.vcf,$(notdir $(basename $(GENOME))))))
PSTOPVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .stop.vcf,$(notdir $(basename $(GENOME))))))
PTFBSVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .tfbs.vcf,$(notdir $(basename $(GENOME))))))

$(POUT)/%.vcf: $(TARGETSDIR)/%.fasta $(REPEATS)	$(GENOME)
	mkdir -p $(POUT)/$(basename $(notdir $<))
	mkdir -p $(POUT)/$(basename $(notdir $<))/input
	cp $(GENOME) $(POUT)/$(basename $(notdir $<))/input/$(notdir $(GENOME))
	cp $< $(POUT)/$(basename $(notdir $<))/input
	-$(PARSNP)/parsnp -r $(POUT)/$(basename $(notdir $<))/input/$(notdir $(GENOME)) -d $(POUT)/$(basename $(notdir $<))/input -p $(PCPU) -v -c -o $(POUT)/$(basename $(notdir $<))/output
	harvesttools -i $(POUT)/$(basename $(notdir $<))/output/parsnp.ggr -V $@.vcf && \
		bedtools subtract -a $@.vcf -b $(REPEATS) > $@.vcf.vcf && \
		$(SRCDIR)/parsnp2vcf $@.vcf.vcf $@ --template $@.vcf && \
		rm $@.vcf && rm $@.vcf.vcf

$(POUT)/%.merged.vcf: $(POUT)/%.vcf $(GENOME)
	cat $< | python2 $(SRCDIR)/merge_variants - $(GENOME) --window 2 > $@	

$(POUT)/%.annotated.vcf: $(POUT)/%.merged.vcf
	$(JAVA7) -jar $(SNPEFFDIR)/snpEff.jar ann $(SNPEFFCHROM) $< > $@

$(POUT)/%.nonsyn.vcf: $(POUT)/%.annotated.vcf
	cat $< | python2 $(SRCDIR)/annvcf2nonsyn - > $@

$(POUT)/%.stop.vcf: $(POUT)/%.annotated.vcf $(GBK)
	cat $< | python2 $(SRCDIR)/annvcf2stops - $(GBK) > $@

$(POUT)/%.tfbs.vcf: $(POUT)/%.vcf $(TFBSTABLE)
	cat $< | python2 $(SRCDIR)/vcf2tfbs $(TFBSTABLE) $(FILESDIR)/pssm - > $@

# Common variants
PCOMMON = $(POUT)/common.tsv
$(PCOMMON): $(PVCFS) $(EVOLUTION)
	$(SRCDIR)/common_variants --exclude $(EVOLUTION) --frequency 0.8 $(PVCFS) > $(PCOMMON)

##############################
## Pairwise reads alignment ##
##############################

# Available reads sets
READS = $(wildcard $(READSDIR)/*)

MVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .vcf,$(notdir $(READ)))))
MANNOTATEDVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .annotated.vcf,$(notdir $(basename $(READ))))))
MSTOPVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .stop.vcf,$(notdir $(basename $(READ))))))
MNONSYNVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .nonsyn.vcf,$(notdir $(READ)))))
MTFBSVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .tfbs.vcf,$(notdir $(READ)))))

GINDEX = $(GENOME).bwt
$(GINDEX): $(GENOME)
	bwa index $(GENOME)

$(MOUT)/%.vcf: $(READSDIR)/% $(GINDEX)
	mkdir -p $(MOUT)/$(basename $(notdir $<))
	mkdir -p $(MOUT)/$(basename $(notdir $<))/subsampled
	mkdir -p $(MOUT)/$(basename $(notdir $<))/aligned
	mkdir -p $(MOUT)/$(basename $(notdir $<))/deduplicated
	mkdir -p $(MOUT)/$(basename $(notdir $<))/realigned
	sample=$$($(SRCDIR)/get_subsample $(GENOME) $$(interleave_pairs $(READSDIR)/$(basename $(notdir $<))/$(READ1) $(READSDIR)/$(basename $(notdir $<))/$(READ2) | count_seqs | awk '{print $$2}') --coverage $(MAXCOVERAGE)) && \
        seqtk sample -s$(SEED) $(READSDIR)/$(basename $(notdir $<))/$(READ1) $$sample > $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ1) && \
	seqtk sample -s$(SEED) $(READSDIR)/$(basename $(notdir $<))/$(READ2) $$sample > $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ2) && \
	bwa mem -t $(MCPU) $(GENOME) $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ1) $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ2) | \
	samtools view -Sb -q 25 -f 2 -F 256 - > $(MOUT)/$(basename $(notdir $<))/aligned/aln.bam && \
	samtools sort $(MOUT)/$(basename $(notdir $<))/aligned/aln.bam $(MOUT)/$(basename $(notdir $<))/aligned/aln.sorted && \
	samtools index $(MOUT)/$(basename $(notdir $<))/aligned/aln.sorted.bam && \
	java -Xmx4g -jar $(PICARDDIR)/MarkDuplicates.jar INPUT=$(MOUT)/$(basename $(notdir $<))/aligned/aln.sorted.bam OUTPUT=$(MOUT)/$(basename $(notdir $<))/deduplicated/aln.dedup.bam METRICS_FILE=$(MOUT)/$(basename $(notdir $<))/deduplicated/metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT && \
	samtools index $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.dedup.bam && \
	java -Xmx$(JAVAMEM)g -jar $(PICARDDIR)/AddOrReplaceReadGroups.jar \
	   I=$(MOUT)/$(basename $(notdir $<))/deduplicated/aln.dedup.bam \
	   O=$(MOUT)/$(basename $(notdir $<))/deduplicated/aln.group.bam \
	   RGPL=illumina RGLB=foo RGPU=run \
	   RGSM=anysample CREATE_INDEX=true && \
	cp $(GENOME) $(MOUT)/$(basename $(notdir $<))/deduplicated/ && \
	java -jar $(PICARDDIR)/CreateSequenceDictionary.jar R=$(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(GENOME)) \
	    O=$(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(basename $(GENOME))).dict \
	    GENOME_ASSEMBLY=genome SPECIES=$(SPECIES) && \
	samtools faidx $(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(GENOME)) && \
	$(JAVA7) -Xmx$(JAVAMEM)g -jar $(GATKDIR)/GenomeAnalysisTK.jar \
	   -T RealignerTargetCreator \
	   -R $(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(GENOME)) \
	   -I $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.group.bam \
	   -o $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.group.bam.intervals && \
	$(JAVA7) -Xmx$(JAVAMEM)g -jar $(GATKDIR)/GenomeAnalysisTK.jar \
	   -T IndelRealigner \
	   -R $(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(GENOME)) \
	   -I $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.group.bam \
	   -targetIntervals $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.group.bam.intervals \
	   -o $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.realn.bam && \
	samtools calmd -brA $(MOUT)/$(basename $(notdir $<))/deduplicated/aln.realn.bam $(MOUT)/$(basename $(notdir $<))/deduplicated/$(notdir $(GENOME)) > $(MOUT)/$(basename $(notdir $<))/realigned/aln.bam && \
	samtools index $(MOUT)/$(basename $(notdir $<))/realigned/aln.bam && \
	freebayes -f $(GENOME) --ploidy $(PLOIDY) --theta $(THETA) --genotype-qualities --standard-filters $(MOUT)/$(basename $(notdir $<))/realigned/aln.bam > $(MOUT)/$(basename $(notdir $<))/raw.vcf && \
	vcffilter $(FILTER) $(MOUT)/$(basename $(notdir $<))/raw.vcf > $@

$(MOUT)/%.annotated.vcf: $(MOUT)/%.vcf
	$(JAVA7) -jar $(SNPEFFDIR)/snpEff.jar ann $(SNPEFFCHROM) $< > $@

$(MOUT)/%.nonsyn.vcf: $(MOUT)/%.annotated.vcf
	cat $< | python2 $(SRCDIR)/annvcf2nonsyn - > $@

$(MOUT)/%.stop.vcf: $(MOUT)/%.annotated.vcf $(GBK)
	cat $< | python2 $(SRCDIR)/annvcf2stops - $(GBK) > $@

$(MOUT)/%.tfbs.vcf: $(MOUT)/%.vcf $(TFBSTABLE)
	cat $< | python2 $(SRCDIR)/vcf2tfbs $(TFBSTABLE) $(FILESDIR)/pssm - > $@

#############################
## Consensus variant calls ##
#############################

CVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(CONSENSUSOUT)/,$(addsuffix .vcf,$(notdir $(basename $(GENOME))))))

$(CONSENSUSOUT)/%.vcf: $(TARGETSDIR)/%.fasta

#######################################
## Gene content variability (LS-BSR) ##
#######################################

$(REFERENCEDIR):
	mkdir -p $(REFERENCEDIR)
$(ALLDIR):
	mkdir -p $(ALLDIR)

# Reference proteome file
REFERENCEFAA = $(REFERENCEDIR)/genome.pep
$(REFERENCEFAA): $(GBK) $(REFERENCEDIR)
	$(SRCDIR)/gbk2faa $(GBK) genome.tmp && \
	$(SRCDIR)/remove_duplicates genome.tmp $(REFERENCEFAA) 

# 1. Reference genes conservation
CONSERVATION = $(REFERENCEDIR)/bsr_matrix_values.txt
$(CONSERVATION): $(REFERENCEFAA) $(REFERENCEDIR)
	mkdir -p $(REFERENCEDIR).tmp && \
		cp $(TARGETSDIR)/*.fasta $(REFERENCEDIR).tmp && \
		cp $(GENOME) $(REFERENCEDIR).tmp
	cd $(REFERENCEDIR) && python2 $(CURDIR)/LS-BSR/ls_bsr.py -d $(REFERENCEDIR).tmp -g $(REFERENCEFAA) -p $(LCPU)
	-rm -rf $(REFERENCEDIR).tmp

# 2. All genes conservation
APPROXPANGENOME = $(ALLDIR)/bsr_matrix_values.txt
$(APPROXPANGENOME): $(REFERENCEFAA) $(ALLDIR)
	mkdir -p $(ALLDIR).tmp
	cp $(TARGETSDIR)/*.fasta $(ALLDIR).tmp
	for genome in $$(cat $(OUTGROUPS)); do rm $(ALLDIR).tmp/$$genome.fasta; done
	cp $(GENOME) $(ALLDIR).tmp
	$(USEARCHDIR)/usearch -cluster_fast $(ALLDIR)/all.faa -id 0.9 -uc $(ALLDIR)/results.uc -centroids $(ALLDIR)/all.pep.tmp && \
	$(SRCDIR)/remove_duplicates $(ALLDIR)/all.pep.tmp $(ALLDIR)/all.pep && \
	cd $(ALLDIR) && python2 $(CURDIR)/LS-BSR/ls_bsr.py -d $(ALLDIR).tmp -g all.pep -p $(LCPU)
	-rm -rf $(ALLDIR).tmp

######################################
## Gene content variability (Roary) ##
######################################

$(RDIR):
	mkdir -p $(RDIR)

$(GFFDIR):
	mkdir -p $(GFFDIR)

REFERENCEPROTEOME = genome.faa
$(REFERENCEPROTEOME): $(GBK)
	$(SRCDIR)/gbk2faa $(GBK) $(REFERENCEPROTEOME) 

REFERENCEGFF = $(GFFDIR)/genome.gff
$(REFERENCEGFF): $(GBK) $(GFFDIR)
	src/gbk2gff $(GBK) $(REFERENCEGFF)

ROARYOUT = $(RDIR)/gene_presence_absence.csv
$(ROARYOUT): $(RDIR) $(REFERENCEGFF) $(GFFS)
	cd $(RDIR) && roary -g 100000 -e -s -v -p $(RCPU) $(GFFDIR)/*.gff
	mkdir -p $(RDIR)/oout
	for strain in $$(cut -d ',' -f 15- $(ROARYOUT) | head -n 1 | sed 's/"//g' | sed 's/,/\t/g' | rev | cut -f 2- | rev); do $(SUBMIT) "$(SRCDIR)/roary2oma $(ROARYOUT) $$strain > $(RDIR)/oout/$$strain.faa.tsv";done

#############################################
## Pairwise gene content variability (OMA) ##
#############################################

# WARNING: this analysis also considers the paralogs

PROTEOMES = $(wildcard $(PROTEOMEDIR)/*)

$(OMAOUT):
	mkdir -p $(OMAOUT)

OTSV = $(foreach PROTEOME,$(PROTEOMES),$(addprefix $(OMAOUT)/,$(addsuffix .tsv,$(notdir $(PROTEOME)))))

OMAPARAMETERS = $(FILESDIR)/parameters.drw

ORTHOXMLLIB = $(SRCDIR)/orthoxml.py
$(ORTHOXMLLIB):
	wget -O $(ORTHOXMLLIB) https://raw.githubusercontent.com/jhcepas/phylogenetic-XML-python-parsers/master/orthoxml.py

$(OMAOUT)/%.tsv: $(PROTEOMEDIR)/% $(REFERENCEFAA) $(ORTHOXMLLIB) $(OMAOUT)
	mkdir -p $(OMAOUT)/$(basename $(notdir $<)) && \
	mkdir -p $(OMAOUT)/$(basename $(notdir $<))/DB && \
	cp $(OMAPARAMETERS) $(OMAOUT)/$(basename $(notdir $<)) && \
	cp $< $(OMAOUT)/$(basename $(notdir $<))/DB/target.fa && \
	cp $(REFERENCEFAA) $(OMAOUT)/$(basename $(notdir $<))/DB/reference.fa && \
	cd $(OMAOUT)/$(basename $(notdir $<)) && $(SUBMIT) "oma -n $(OCPU) && \
       	python2 $(SRCDIR)/omah2tsv $(OMAOUT)/$(basename $(notdir $<))/OutputFolder/HierarchicalGroups.orthoxml $@"

####################
## K-mer counting ##
####################

GENOMES = $(wildcard $(TARGETSDIR)/*.fasta)
NAMES = $(foreach GENOME,$(GENOMES),$(notdir $(basename $(GENOME))))
KMERS = $(foreach GENOME,$(GENOMES),$(addprefix $(KMOUT)/mers/,$(addsuffix .jf,$(notdir $(basename $(GENOME))))))
KCOUNTS = $(foreach GENOME,$(GENOMES),$(addprefix $(KMOUT)/counts/,$(addsuffix .txt,$(notdir $(basename $(GENOME))))))

$(KMOUT):
	mkdir -p $(KMOUT) && \
	mkdir -p $(KMOUT)/mers && \
	mkdir -p $(KMOUT)/counts

$(KMOUT)/mers/%.jf: $(TARGETSDIR)/%.fasta $(KMOUT)
	$(JELLYFISHDIR)/jellyfish count -C -m $(WORD) -s $(KHASH)M $< -o $@

REFERENCEKMER = $(KMOUT)/genome.jf
$(REFERENCEKMER): $(KMOUT)
	$(JELLYFISHDIR)/jellyfish count -C -m $(WORD) -s $(KHASH)M $(GENOME) -o $(REFERENCEKMER)

ALLGENOMES = all.fasta
$(ALLGENOMES): $(GENOME)
	cat $(GENOME) $(TARGETSDIR)/*.fasta > $(ALLGENOMES)

$(KMOUT)/counts/%.txt: $(TARGETSDIR)/%.fasta $(KMOUT) $(ALLGENOMES) $(KMERS)
	$(JELLYFISHDIR)/jellyfish query $(KMOUT)/mers/$(basename $(notdir $<)).jf -s $(ALLGENOMES) | cut -d" " -f2 > $@

REFERENCECOUNT = $(KMOUT)/genome.txt
$(REFERENCECOUNT): $(KMOUT) $(ALLGENOMES) $(REFERENCEKMER)
	$(JELLYFISHDIR)/jellyfish query $(REFERENCEKMER) -s $(ALLGENOMES) > $(REFERENCECOUNT)

KTMP = $(KMOUT)/kmers.txt
KTABLE = $(KTMP).gz
$(KTABLE): $(KCOUNTS) $(REFERENCECOUNT)
	echo -e "word" $(notdir $(basename $(REFERENCECOUNT))) $(NAMES) > $(KTMP) && \
	paste $(REFERENCECOUNT) $(KCOUNTS) >> $(KTMP) && \
	gzip $(KTMP) && \
	zcat $(KTABLE) | tail -n+2 | split -l 10000000 - $(KMOUT)/kmer_split_ && \
	for i in $(ls $(KMOUT)/kmer_split_*); do sed -i "1i$$(zcat $(KTABLE) | head -n 1)" $i; done && \
	gzip $(KMOUT)/kmer_split_*

###############################
## K-mer counting (fsm-lite) ##
###############################

$(KMLITEOUT):
	mkdir -p $(KMLITEOUT)

KLITETABLE = $(KMLITEOUT)/kmers.txt.gz
$(KLITETABLE): $(KMLITEOUT)
	-rm $(KMLITEOUT)/input.txt 
	for infile in $$(ls $(TARGETSDIR) | grep ".fasta" | awk -F'_' '{print $$1}' | sort | uniq); do \
	  echo -e $$infile"\t"$(TARGETSDIR)/$$(ls $(TARGETSDIR) | grep $$infile | sort | uniq | tail -n 1) >> $(KMLITEOUT)/input.txt; \
	done
	fsm-lite -l $(KMLITEOUT)/input.txt -t $(KMLITEOUT)/tmp.txt -m 9 -M 100 -f 1 -s 1 -S 99 -v > $(KMLITEOUT)/kmers.txt && \
	split -d -n l/16 $(KMLITEOUT)/kmers.txt $(KMLITEOUT)/fsm_out && \
	gzip $(KMLITEOUT)/fsm_out* && \
	gzip $(KMLITEOUT)/kmers.txt

#################################
## Strains tree (using parsnp) ##
#################################

$(TOUT):
	mkdir -p $(TOUT)

TREE = $(TOUT)/output/parsnp.tree
$(TREE): $(GENOME) $(TOUT)
	mkdir -p $(TOUT)/input && \
	cp $(GENOME) $(TOUT)/input/ && \
	cp $(TARGETSDIR)/*.fasta $(TOUT)/input && \
	$(PARSNP)/parsnp -r $(TOUT)/input/$(notdir $(GENOME)) -d $(TOUT)/input -p $(PCPU) -v -c -o $(TOUT)/output && \
	rm -rf $(TOUT)/input/

###################################
## Strains distance (using mash) ##
###################################

$(MASHDIR):
	mkdir -p $(MASHDIR)

MASHDISTANCES = $(MASHDIR)/distances.tab
$(MASHDISTANCES): $(GENOME) $(MASHDIR)
	cd $(MASHDIR) && \
	echo ../$(notdir $(GENOME)) > genomes.txt && \
	ls $(TARGETSDIR)/*.fasta >> genomes.txt && \
	mash sketch -l genomes.txt -o genomes -p $(MCPU) && \
	mash dist genomes.msh genomes.msh -p $(MCPU) | python -c "import os;import sys;print('\n'.join([os.path.split(l.split('\t')[0])[-1].split('_')[0].split('.')[0] + '\t' + os.path.split(l.split('\t')[1])[-1].split('_')[0].split('.')[0] + '\t' + '\t'.join(l.rstrip().split('\t')[2:]) for l in sys.stdin]))" > $(MASHDISTANCES)

####################
## RegulonDB data ##
####################

PSSMDIR = $(FILESDIR)/pssm
$(PSSMDIR):
	mkdir -p $(PSSMDIR)

TFBS = $(FILESDIR)/BindingSiteSet.txt
$(TFBS):
	wget -O $(TFBS) http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt

PSSM = $(FILESDIR)/PSSMSet.txt
$(PSSM): $(PSSMDIR)
	wget -O $(PSSM) http://regulondb.ccg.unam.mx/menu/download/datasets/files/PSSMSet.txt
	$(SRCDIR)/retrieve_pssm $(PSSM) $(PSSMDIR)

TFBSTABLE = $(FILESDIR)/tfbs.txt
$(TFBSTABLE): $(TFBS) $(PSSM) $(GBK)
	$(SRCDIR)/get_regulated_genes $(GBK) $(TFBS) --details > $(TFBSTABLE).tmp
	$(SRCDIR)/correct_tfbstable $(TFBSTABLE).tmp $(PSSMDIR) $(GBK) > $(TFBSTABLE)

########################
## gkm-SVM generation ##
########################

POSITIVESET = $(GKDIR)/positive.fasta
NEGATIVESET = $(GKDIR)/negative.fasta

$(POSITIVESET): $(GBK) $(TFBSTABLE)
	$(SRCDIR)/generate_sets $(GBK) $(TFBSTABLE) 40 $(POSITIVESET) $(NEGATIVESET) --intergenic

KERNEL = $(GKDIR)/kernel.out
$(KERNEL): $(POSITIVESET)
	$(GKDIR)/gkmsvm_kernel -d 3 $(POSITIVESET) $(NEGATIVESET) $(KERNEL)

SVMALPHA = $(GKDIR)/svmtrain_svalpha.out
SVSEQ = $(GKDIR)/svmtrain_svseq.fa
$(SVMALPHA): $(KERNEL) $(POSITIVESET) 
	$(GKDIR)/gkmsvm_train $(KERNEL) $(POSITIVESET) $(NEGATIVESET) $(GKDIR)/$(shell basename $(SVMALPHA) _svalpha.out)

WEIGHTS = $(GKDIR)/weights.txt
$(WEIGHTS): $(SVMALPHA)
	$(SRCDIR)/nrkmers 10 $(GKDIR)/nr10.fasta
	$(GKDIR)/gkmsvm_classify -d 3 $(GKDIR)/nr10.fasta $(SVSEQ) $(SVMALPHA) $(GKDIR)/nr10_scores.txt
	$(SRCDIR)/generate_weights $(GKDIR)/nr10_scores.txt > $(WEIGHTS)

#############
## Targets ##
#############

all: ksnp parsnp map conservation oma kmers roary
ksnp: $(KVCFS)
parsnp: $(PVCFS) $(PANNOTATEDVCFS) $(PMERGEDVCFS)
map: $(MVCFS) $(MANNOTATEDVCFS)
common: $(PCOMMON)
consensus: $(CVCFS) $(MVCFS) $(PVCFS) $(KVCFS)
conservation: $(CONSERVATION) $(APPROXPANGENOME)
oma: $(OTSV)
kmers: $(KTABLE)
fsm: $(KLITETABLE)
roary: $(ROARYOUT)
tree: $(TREE) $(MASHDISTANCES)
nonsyn: $(PNONSYNVCFS)
tfbs: $(MTFBSVCFS) $(PTFBSVCFS)
stop: $(PSTOPVCFS) 
regulondb: $(TFBSTABLE)
gksvm: $(WEIGHTS)
pangenome: $(RPANGENOME)
snps: $(RSNPSM)
stops: $(RSTOP)

.PHONY: all ksnp parsnp map common conservation oma kmers fsm roary tree nonsyn tfbs stop regulondb gksvm pangenome snps stops
