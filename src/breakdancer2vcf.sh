#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
for i in 01_SV_calling/*
do
id=$(basename $i)

python ~/software/script/WGS_SV/src/breakdancer2vcf.py -i $i/breakdancer/$id.cfg.SV.output -o $i/breakdancer/$id.cfg.SV.vcf    
echo "$i completed"
done

