ecoliNTdb/src/get_strains_genomes --db ecoliNTdb/ntdb.sqlite | src/get_genomes_id > genomes.start.txt
mkdir -p genomes
mkdir -p proteomes
mkdir -p gff
mkdir -p prokka
for genome in $(awk '{print $1}' genomes.start.txt)
do
    echo $genome;
    src/gbk2fasta ecoliNTdb/genomes/$(grep $genome genomes.start.txt | awk '{print $2}') genomes/$genome.fasta;
done

for i in $(ls genomes);
do
  strain=$(basename $i .fasta);
  prokka --outdir prokka/$strain --force --prefix $strain --addgenes --locustag $strain --mincontiglen 200 --genus Escherichia -species coli --strain $strain --norrna --notrna --proteins genome.faa genomes/$i;
done

cp prokka/*/*.faa proteomes/
cp prokka/*/*.gff gff/

for i in $(cat outgroups.txt);
do
  echo $i;
  rm genomes/$i*;
  rm gff/$i*;
  rm proteomes/$i*;
done
