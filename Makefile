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
# Directory where the parsnp binary is
PARSNP = Parsnp-Linux64-v1.2

# Parameters
# kSNP CPUs
KCPU = 20
# NOTE: optimal k-mer shoudld be derived from
# a Kchooser run
KKMER = 19
# parsnp CPUs
PCPU = 10

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

VCFS = $(foreach GENOME,$(GENOMES),$(addprefix $(POUT)/,$(addsuffix .vcf,$(notdir $(GENOME)))))

MASKEDGENOME = genome.masked.fasta
$(MASKEDGENOME): $(GENOME) $(REPEATS)
	bedtools maskfasta -fi $(GENOME) -bed $(REPEATS) -fo $(MASKEDGENOME)

$(POUT)/%.vcf: $(TARGETSDIR)/% $(MASKEDGENOME)	
	mkdir -p $(POUT)/$(basename $(notdir $<))
	mkdir -p $(POUT)/$(basename $(notdir $<))/input
	cp $(MASKEDGENOME) $(POUT)/$(basename $(notdir $<))/input
	cp $< $(POUT)/$(basename $(notdir $<))/input
	$(PARSNP)/parsnp -r genomes/$(GENOME) -d genomes -p $(PCPU) -v -c -o $(POUT)/$(basename $(notdir $<))/output
	harvesttools -i $(POUT)/$(basename $(notdir $<))/output/parsnp.ggr -V $@

all: ksnp
ksnp: $(KOUTPUT)
parsnp: $(VCFS)

.PHONY: all ksnp parsnp
