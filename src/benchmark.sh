#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
for vcf in *_DEL.vcf.gz
do
id=$(basename $vcf .vcf.gz)
# vcf=pbsv_DEL.vcf.gz

baseline=/WORK2/gzfezx_shhli_3/BioDatahub/GIAB/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6_DEL.vcf.gz
base_bed=/WORK2/gzfezx_shhli_3/BioDatahub/GIAB/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
truvari=/BIGDATA2/gzfezx_shhli_2/software/bin/truvari
$truvari bench -b $baseline -c $vcf --includebed $base_bed  --pctseq 0 -o $id -r 1000
done


for vcf in *_INS.vcf.gz
do
id=$(basename $vcf .vcf.gz)
# vcf=pbsv_DEL.vcf.gz

baseline=/WORK2/gzfezx_shhli_3/BioDatahub/GIAB/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6_INS.vcf.gz
base_bed=/WORK2/gzfezx_shhli_3/BioDatahub/GIAB/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
truvari=/BIGDATA2/gzfezx_shhli_2/software/bin/truvari
$truvari bench -b $baseline -c $vcf --includebed $base_bed  --pctseq 0 -o $id -r 1000
done
