# Input files
# The genome file should be put in the same directory as the makefile
GENOME = genome.fasta
GBK = genome.gbk
# Directory containing the target genomes fasta files
TARGETSDIR = $(CURDIR)/genomes
# Output directory for kSNP
KOUT = $(CURDIR)/kout

# Directories and parameters
# kSNP CPUS
KCPU = 10
# NOTE: optimal k-mer shoudld be derived from
# a Kchooser run
KKMER = 19

# Anything below this point should not be changed
$(TARGETSDIR): $(GENOME)
	mkdir -p $(TARGETSDIR)

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

all: ksnp
ksnp: $(KOUTPUT)

.PHONY: all ksnp
