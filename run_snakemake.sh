#!/bin/bash
# yhbatch -p e9 -n 1
# yhbatch -N 1 -n 64 -p deimos
source /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/activate manta2
pwd=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240801_89_SVs/test/15100418BFF2A
smk=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/run.all.smk
threads=64
/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake -d $pwd -c $threads -pk --configfile ${pwd}/config.yaml --default-resources tmpdir="\"${pwd}\"" -s ${smk} --rerun-triggers mtime > ${pwd}/snakemake.err.txt 2>&1 

# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)" xx.email
