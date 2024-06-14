#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 
echo "##################   start run delly   #######################"

mkdir delly
cd delly
bam=../00116021276M31BFF2A.sorted.markdup.BQSR.cram
sample=$(basename $bam .sorted.markdup.BQSR.cram)
ref=/WORK/gzfezx_shhli_3/BioDatahub/human_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
delly call -o delly.bcf -g $ref $bam

bcftools view delly.bcf > delly.vcf

# delly cnv -g example/ref.fa -m example/map.fa.gz -c out.cov.gz -o cnv.bcf example/sr.bam
echo "##################   end run delly   #######################"
