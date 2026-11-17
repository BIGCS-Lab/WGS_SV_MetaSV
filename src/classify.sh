
for i in breakdancer.vcf.gz cnvnator.vcf.gz  lumpy.vcf.gz  manta.vcf.gz  metasv.vcf.gz
do
# i=breakdancer.vcf.gz 
id=$(basename $i .vcf.gz)
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/python \
/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/split_vcf_bysvtype.py $i $id

if [ -e ${id}_DUP.vcf ] && [ -e ${id}_INS.vcf ] ; then
    /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/python \
    /BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/cat_DUP_INS.py \
        -i ${id}_INS.vcf -d ${id}_DUP.vcf -o ${id}_DUP_INS.vcf
else
    if [ -e ${id}_INS.vcf ]; then
    cp -if  ${id}_INS.vcf  ${id}_DUP_INS.vcf
    else
    cp -if  ${id}_DUP.vcf  ${id}_DUP_INS.vcf
    fi
fi

/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/python \
/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/change_vcf_info.py \
${id}_DEL.vcf tmp1 && mv -if tmp1 ${id}_DEL.vcf
bgzip -c ${id}_DEL.vcf > ${id}_DEL.vcf.gz
tabix ${id}_DEL.vcf.gz


cat /BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/vcff_hg37.header ${id}_DUP_INS.vcf > tmp && mv -if tmp ${id}_DUP_INS.vcf

/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/python \
/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/change_vcf_info.py \
${id}_DUP_INS.vcf tmp1 && mv -if tmp1 ${id}_DUP_INS.vcf

/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/python \
/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/src/change_vcf_info.py \
${id}_DEL.vcf tmp1 && mv -if tmp1 ${id}_DEL.vcf

bcftools sort ${id}_DUP_INS.vcf|bgzip -c > ${id}_DUP_INS.vcf.gz
tabix ${id}_DUP_INS.vcf.gz


done
