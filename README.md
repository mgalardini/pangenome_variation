pangenome_variation
===================

Scripts and pipeline to inspect genetic variation in a series of bacterial genomes

Note
----

The pipeline and scripts come with limited documentation.
Please do get in touch with the author (Marco Galardini, marco@ebi.ac.uk) if you need any guidance.

Usage
-----

A reference genome in FASTA and Genbank format is needed
(deafult filenames are genome.fasta and genome.gbk).
All the genomes to be analysed should be assemblies: place 
nucleotides fasta files in the `genomes` directory
(`genomes/*.fasta`), protein fasta files in the `proteomes` directory
(`proteomes/*.faa`) and gff files in the `gff` directory (`gff/*.gff`).
We reccommend using prokka to generate the `.faa` and `.gff` files.

The makefile contains the various bits of the pipeline:
* `make tree`: core genome alignment phylogenetic tree and mash whole genome kmer distance
* `make roary`: pangenome
* `make oma`: pairwise pangenome for each strain against the reference
* `make nonsyn stop`: pairwise alignment of each strain agains the reference to derive SNPs

You might want to type `make -n TARGET` first to make sure which commands are gonna be launched

(minimum) prerequisites
-------------

* prokka
* parsnp and harvest
* mash
* snpeff
* roary
* oma
* python (2.7+ AND 3.3+), plus the following libraries:
    * biopython
    * bcbio-gff
    * numpy
    * pandas
    * pyvcf

Copyright
---------

Copyright (C) <2015> EMBL-European Bioinformatics Institute

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
GNU General Public License for more details.

Neither the institution name nor the name pangenome_variation
can be used to endorse or promote products derived from
this software without prior written permission.
For written permission, please contact <marco@ebi.ac.uk>.

Products derived from this software may not be called pangenome_variation
nor may pangenome_variation appear in their names without prior written
permission of the developers. You should have received a copy
of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.
