#!/usr/bin/env bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta
smk=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV_mul/run.all.smk
id="w1"

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`

/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake \
        --resources mem_mb=120000 --configfile $id.config.yaml \
        -j 24 -pk -s ${smk} 2> $id.snakemake.err.txt

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
sumTime=$[ $endTime_s - $startTime_s ]
echo "$startTime ---> $endTime" "Total:$sumTime seconds"
