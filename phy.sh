#!/bin/bash
mkdir vcf/merge
vcf-merge vcf/fin_vcf/* > vcf/merge/merged.vcf
python vcf2phylip.py --input merged.vcf
python ascbias.py -p merged.min4.phy --outfile merged_out.phy
cat merged_out.phy | sed -e "s/*/-/g" > merged_out_2.phy
modeltest-ng -i merged_out_2.phy
#raxml-ng --msa merged_out_2.phy --model GTR --all --bs-trees 10000 --threads 16
raxml-ng --msa merged_out_2.phy --model GTR+ASC_LEWIS --all --bs-trees 10000 --threads 25
