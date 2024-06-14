#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta
env="/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta/bin"
bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
REF="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
WD=`pwd`
readlength=100
lumpyinput=../lumpy/00114031182M22BFF2.lumpy.vcf
breakdancerinput=../breakdancer/00114031182M22BFF2.cfg.SV.output
cnvnatorinput=../cnvnator/00114031182M22BFF2.vcf
mantainput=../manta/results/variants/diploidSV.vcf.gz
${env}/python ${env}/run_metasv.py \
            --reference ${REF} --sample ${sample} --disable_assembly --num_threads 24 \
            --enable_per_tool_output --keep_standard_contigs --mean_read_length ${readlength} \
            --outdir ./ --workdir ${WD}/tmp_work \
            --breakdancer_native ${breakdancerinput} \
            --manta_vcf ${mantainput} \
            --cnvnator_vcf ${cnvnatorinput} \
            --lumpy_vcf ${lumpyinput} &>${WD}/${sample}.metasv.log
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta

bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/breakseq2/bin/perl /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/breakseq2/bin/bam2cfg.pl -g -h $bam >$sample.cfg
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/breakseq2/bin/breakdancer-max $sample.cfg > $sample.cfg.SV.output 
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/breakseq2

bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
bam2cfg.pl -g -h $bam >$sample.cfg
breakdancer-max $sample.cfg > $sample.cfg.SV.output 
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta

bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
REF="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

pindel -i pindel.conf -f $REF -o $sample
# pindel2vcf -r <参考基因组文件> -R <参考基因组名称> -d 参考基因组日期 -p <pindel输出文件> -e <最小的reads数>
pindel2vcf -r $REF -R hg38 -d 20140111 -P $sample -v $sample.vcf
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta
env="/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta/bin"
bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
REF="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
WD=`pwd`

python ${env}/configManta.py --bam ${bam} --referenceFasta ${REF} --runDir ${WD}/ &>/dev/null
python ${WD}/runWorkflow.py -j 24 &>/dev/null
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta

bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
REF="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
samtools view -@ 24 -b -F 1294 ${bam} > ${sample}.discordants.unsorted.bam 
samtools view -h ${bam} | extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${sample}.splitters.unsorted.bam
samtools sort -@ 24 ${sample}.discordants.unsorted.bam -o ${sample}.discordants.bam
samtools sort -@ 24 ${sample}.splitters.unsorted.bam -o ${sample}.splitters.bam
lumpyexpress -B ${bam} -S ${sample}.splitters.bam -D ${sample}.discordants.bam -o ${sample}.lumpy.vcf > ${sample}.lumpy.runlog
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator

bam="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/00114031182M22BFF2/00114031182M22BFF2.bam"
sample=$(basename $bam .bam)
REF="/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
# bin_size 根据测序深度的不同进行调整：100X: 30；20-30X: 100；4-6X: 500；再低就1000。
bin_size=500

# Step1  从bam文件中提取比对信息
cnvnator -root ${sample}.root -tree $bam -chrom $(seq -f 'chr%g' 1 22) chrX chrY
# Step2生成柱形图 
cnvnator -root ${sample}.root -his ${bin_size} -fasta $REF
# Step3统计量计算
cnvnator -root ${sample}.root -stat ${bin_size} 
Step4 RD信号分割
cnvnator -root ${sample}.root -partition ${bin_size} -ngc
# Step5 call CNV
cnvnator -root ${sample}.root -call ${bin_size} -ngc > ${sample}.cnvnator
# 把sample.cnvnator文件 转换成vcf文件
cnvnator2VCF.pl -prefix ${sample} -reference $REF \
    ${sample}.cnvnator /BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/GRCh38/chr > ${sample}.vcf
#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate manta
cram=00114031182M22BFF2.sorted.markdup.BQSR.cram
sample=$(basename $cram .sorted.markdup.BQSR.cram)
samtools view -@ 24 -bh $cram -o ${sample}.bam
samtools index -@ 24 ${sample}.bam
