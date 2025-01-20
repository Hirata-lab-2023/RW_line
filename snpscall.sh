#!/bin/bash
for t in `cat direc2.txt`
do
mkdir ${t}
done
for t in `cat SRA.txt`
do
fastq-dump --gzip --split-files --outdir rawdata ${t} &
done
wait
for t in `cat SRA.txt`
do
i=($(grep -n "${t}" SRA.txt | awk -F ':' '{print $1}'))
# trimming
fastp --thread 8 --in1 rawdata/${t}_1.fastq.gz --in2 rawdata/${t}_2.fastq.gz --out1 trimed/0${i}_1_trim.fastq.gz --out2 trimed/0${i}_2_trim.fastq.gz --detect_adapter_for_pe &
done
wait
for i in {1..36}; do bwa mem -t 5 Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa trimed/0${i}_1_trim.fastq.gz  trimed/0${i}_2_trim.fastq.gz -o sam/0${i}.sam & done
wait
rm trimed/*
for i in {1..36}
do
rm trimed/0${i}_*_trim.fastq.gz &
#
# sam > bam
samtools view -h sam/0${i}.sam | awk '$17 !~ /XA:/|| $1 ~ /^@/' | samtools view -bS - > bam/0${i}.uniq.bam &
done
wait
for i in {1..36}
do
rm sam/0${i}.sam
samtools sort bam/0${i}.uniq.bam > bam/0${i}.uniq.sort.bam &
done
wait
for i in {1..36}
do
java -Xmx50G -jar picard.jar AddOrReplaceReadGroups I=bam/0${i}.uniq.sort.bam O=bam/0${i}.uniq.sort.addG.bam RGID=0${i}  RGLB=0${i} RGPL=BGI RGPU=run_barcode RGSM=0${i} &
done
wait
for i in {1..36}
do
rm bam/0${i}.uniq.sort.bam
samtools index bam/0${i}.uniq.sort.addG.bam &
done
wait
for i in {1..36}
do
gatk MarkDuplicates -I bam/0${i}.uniq.sort.addG.bam -O bam/0${i}.uniq.sort.addG_markdup.bam -M bam/0${i}.uniq.sort.addG_markdup_matrics.txt
samtools index bam/0${i}.uniq.sort.addG_markdup.bam &
done
wait
for i in {1..36}
do
rm bam/0${i}.uniq.sort.addG.bam
done
# bam > vcf
for i in {1..36}
do
mkdir vcf/0${i}
for t in `cat chr_mt.txt`
do
gatk --java-options "-Xmx4g" HaplotypeCaller -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa -A StrandBiasBySample -mbq 20 -L ${t} -I bam/0${i}.uniq.sort.addG_markdup.bam -O vcf/0${i}/0${i}_${t}.vcf &
done
done
wait
for i in {1..36}
do
for t in `cat chr_mt.txt`
do
bgzip vcf/0${i}/0${i}_${t}.vcf &
done
wait
for t in `cat chr_mt.txt`
do
tabix -p vcf vcf/0${i}/0${i}_${t}.vcf.gz &
done
done
wait
for i in {1..36}
do
for t in `cat chr_mt.txt`
do
gatk --java-options "-Xmx4g" SelectVariants -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa -V vcf/0${i}/0${i}_${t}.vcf.gz --select-type-to-include SNP -O vcf/0${i}/0${i}_${t}.snps.vcf &
done
done
wait
for i in {1..36}
do
for t in `cat chr_mt.txt`
do
gatk --java-options "-Xmx4g" VariantFiltration -R Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa -V vcf/0${i}/0${i}_${t}.snps.vcf -O vcf/0${i}/0${i}_${t}.snps.MQ60.vcf -filter "QD < 2.0" --filter-name "QD_2.0" -filter "FS > 60.0" --filter-name "FS_60" -filter "MQ < 60.0" --filter-name "MQ_60" -filter "MQRankSum < -12.5" --filter-name "MQRankSum_-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_-8.0" -G-filter "DP < 4" --G-filter-name "lowiDP" &
done
done
wait
for i in {1..36}
do
less vcf/0${i}/0${i}_1.snps.MQ60.vcf | grep -e "#" > vcf/0${i}/0${i}_1.head.snps.MQ60.vcf
for t in `cat chr.txt`
do
less vcf/0${i}/0${i}_${t}.snps.MQ60.vcf | grep -v "#" | grep -e "PASS" > vcf/0${i}/0${i}_${t}_cut.snps.MQ60.vcf
done
wait
bgzip -f vcf/0${i}/0${i}_MT_cut.snps.MQ60.vcf
cat vcf/0${i}/0${i}_1.head.snps.MQ60.vcf vcf/0${i}/0${i}_*_cut.snps.MQ60.vcf > vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf
bgzip -f vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf
vcf-sort vcf/0${i}/0${i}_all_cat_snps.MQ60.vcf.gz | bgzip -c > vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz
tabix -p vcf vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz
java -jar snpEff.jar -v -stats vcf/0${i}/0${i}_ano_all.html GRCz11_109 vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz > vcf/0${i}/0${i}_ano_all.vcf
rm vcf/0${i}/0${i}_*_cut.snps.MQ60.vcf
rm vcf/0${i}/0${i}_1.head.snps.MQ60.vcf
cat vcf/0${i}/0${i}_ano_all.vcf | java -jar SnpSift.jar filter "( ANN[*].IMPACT == 'HIGH')" > vcf/0${i}/0${i}_ano_all_HIGH.vcf &
cat vcf/0${i}/0${i}_ano_all.vcf | java -jar SnpSift.jar filter "( ANN[*].IMPACT == 'MODERATE')" > vcf/0${i}/0${i}_ano_all_MOD.vcf &
done
wait
mkdir vcf/HIGH
mkdir vcf/MOD
for i in {1..36}
do
bgzip -f vcf/0${i}/0${i}_ano_all.vcf
tabix -p vcf vcf/0${i}/0${i}_ano_all.vcf.gz&
bgzip -f vcf/0${i}/0${i}_ano_all_HIGH.vcf
tabix -p vcf vcf/0${i}/0${i}_ano_all_HIGH.vcf.gz
mv vcf/0${i}/0${i}_ano_all_HIGH.vcf* vcf/HIGH&
bgzip -f vcf/0${i}/0${i}_ano_all_MOD.vcf
tabix -p vcf vcf/0${i}/0${i}_ano_all_MOD.vcf.gz
mv vcf/0${i}/0${i}_ano_all_MOD.vcf* vcf/MOD&
done
for i in {1..36}
do
mkdir vcf/fin_vcf
mv vcf/0${i}/0${i}_all_sort_cat_snps.MQ60.vcf.gz vcf/fin_vcf/
done
