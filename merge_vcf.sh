#!/bin/bash

mkdir all_vcf
find pout -maxdepth 1 -type f -name '*.vcf' ! -name '*.merged.*' ! -name '*.nonsyn.*' ! -name '*.stop.*' ! -name '*.annotated.*' -exec cp {} all_vcf/ \;
cd all_vcf
for i in $(find . -name '*.vcf'); do bgzip $i && tabix -p vcf $i.gz; done
vcf-merge *.vcf.gz > all.vcf
rm keep.txt
for i in $(find . -type f -name '*.vcf.gz'); do echo $(echo $i | sed 's/.\///g' | sed 's/.vcf.gz//g').fasta >> keep.txt; done
cat all.vcf | vcf-subset -r -e -c keep.txt -f > subset.vcf
../src/vcf_filter subset.vcf | grep -v "##source" | grep -v "##INFO" > snps.vcf
bgzip snps.vcf
rm subset.vcf
rm all.vcf
rm NT*vcf*
