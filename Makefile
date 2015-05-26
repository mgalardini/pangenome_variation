# Input files
# The genome file should be put in the same directory as the makefile
GENOME = genome.fasta
GBK = genome.gbk
# Directory containing the target genomes fasta files
TARGETSDIR = $(CURDIR)/genomes
# Output directory for kSNP
KOUT = $(CURDIR)/kout
# Output directory for parsnp
POUT = $(CURDIR)/pout
# Output directory for reads alignment
MOUT = $(CURDIR)/mout
# Directory where the parsnp binary is
PARSNP = $(SOFTDIR)/Parsnp-Linux64-v1.2
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
READ1 = READ1.tar.gz
READ2 = READ2.tar.gz

# Parameters
# kSNP CPUs
KCPU = 20
# NOTE: optimal k-mer shoudld be derived from
# a Kchooser run
KKMER = 19
# parsnp CPUs
PCPU = 1
# Reads mapping CPUs
MCPU = 20
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

# Anything below this point should not be changed
$(TARGETSDIR): $(GENOME)
	mkdir -p $(TARGETSDIR)

# kSNP
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

# Pairwise parsnp
GENOMES = $(wildcard $(TARGETSDIR)/*)

$(POUT):
	mkdir -p $(POUT)

# Mask repeats in the reference genome
REPEATS = repeats.bed
$(REPEATS): $(GENOME)
	nucmer --maxmatch --nosimplify $(GENOME) $(GENOME) && \
	show-coords -r -T out.delta -H | tail -n+2 > repeats.txt && \
	awk '{print $$8"\t"$$1"\t"$$2}' repeats.txt > $(REPEATS)

PVCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .vcf,$(notdir $(GENOME)))))

MASKEDGENOME = genome.masked.fasta
$(MASKEDGENOME): $(GENOME) $(REPEATS)
	bedtools maskfasta -fi $(GENOME) -bed $(REPEATS) -fo $(MASKEDGENOME)

$(POUT)/%.vcf: $(TARGETSDIR)/% $(MASKEDGENOME)	
	mkdir -p $(POUT)/$(basename $(notdir $<))
	mkdir -p $(POUT)/$(basename $(notdir $<))/input
	cp $(MASKEDGENOME) $(POUT)/$(basename $(notdir $<))/input/$(notdir $(GENOME))
	cp $< $(POUT)/$(basename $(notdir $<))/input
	$(PARSNP)/parsnp -r $(POUT)/$(basename $(notdir $<))/input/$(notdir $(GENOME)) -d $(POUT)/$(basename $(notdir $<))/input -p $(PCPU) -v -c -o $(POUT)/$(basename $(notdir $<))/output
	harvesttools -i $(POUT)/$(basename $(notdir $<))/output/parsnp.ggr -V $@

# Pairwise reads alignment

# Available reads sets
READS = $(sort $(dir $(wildcard $(READSDIR)/*/)))

MVCFS = $(foreach READ,$(READS),$(addprefix $(MOUT)/,$(addsuffix .vcf,$(notdir $(READ)))))

GINDEX = $(GENOME).bwt
$(GINDEX): $(GENOME)
	bwa index $(GENOME)

$(MOUT)/%.vcf: $(READSDIR)/%
	mkdir -p $(MOUT)/$(basename $(notdir $<))
	mkdir -p $(MOUT)/$(basename $(notdir $<))/trimmed
	mkdir -p $(MOUT)/$(basename $(notdir $<))/subsampled
	mkdir -p $(MOUT)/$(basename $(notdir $<))/aligned
	mkdir -p $(MOUT)/$(basename $(notdir $<))/deduplicated
	mkdir -p $(MOUT)/$(basename $(notdir $<))/realigned
	interleave_pairs $(READSDIR)/$(basename $(notdir $<))/$(READ1) $(READSDIR)/$(basename $(notdir $<))/$(READ2) | \
	trim_edges -l 9 --paired_reads | \
	deinterleave_pairs -z -o $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ1) $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ2) && \
	sample=$$($(SRCDIR)/get_subsample $(GENOME) $$(interleave_pairs $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ1) $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ2) | count_seqs | awk '{print $$2}') --coverage $(MAXCOVERAGE)) && \
        seqtk sample -s$(SEED) $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ1) $$sample > $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ1) && \
	seqtk sample -s$(SEED) $(MOUT)/$(basename $(notdir $<))/trimmed/$(READ2) $$sample > $(MOUT)/$(basename $(notdir $<))/subsampled/$(READ2) && \
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
	    GENOME_ASSEMBLY=genome SPECIES=$(SPECIES)
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
	
all: ksnp parsnp map
ksnp: $(KOUTPUT)
parsnp: $(PVCFS)
map: $(MVCFS) 

.PHONY: all ksnp parsnp map
