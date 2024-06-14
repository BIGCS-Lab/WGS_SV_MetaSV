#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 
mkdir manta && cd manta
bam=00116021276M31BFF2A.sorted.markdup.BQSR.cram
sample=$(basename $bam .sorted.markdup.BQSR.cram)
ref=/WORK/gzfezx_shhli_3/BioDatahub/human_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
configManta.py --bam $bam \
--referenceFasta  $ref --runDir manta
manta/runWorkflow.py
