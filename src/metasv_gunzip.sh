#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
for i in 01_SV_calling/*
do
id=$(basename $i)

gunzip -c $i/metasv/$id.SV.vcf.gz > $i/metasv/$id.SV.vcf

echo "$i completed"
done
