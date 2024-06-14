#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 
mkdir gridss && cd gridss
bam=../00116021276M31BFF2A.sorted.markdup.BQSR.cram
sample=$(basename $bam .sorted.markdup.BQSR.cram)
ref=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
gridss \
  -r $ref \
  -o gridss.vcf \
  -t 24 \
  $bam
