$(REPODIR)ecoliNTdb/src/get_strains_genomes --db $(REPODIR)ecoliNTdb/ntdb.sqlite | src/get_genomes_id > genomes.start.txt
mkdir -p genomes
mkdir -p proteomes
for genome in $(awk '{print $1}' genomes.start.txt)
do
    echo $genome;
    src/gbk2fasta $(REPODIR)ecoliNTdb/genomes/$(grep $genome genomes.start.txt | awk '{print $2}') genomes/$genome.fasta;
    src/gbk2faa $(REPODIR)ecoliNTdb/genomes/$(grep $genome genomes.start.txt | awk '{print $2}') proteomes/$genome.faa;
done
