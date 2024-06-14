#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 
echo "##################   start run lumpy   #######################"
mkdir lumpy
cd lumpy
bam=../00116021276M31BFF2A.sorted.markdup.BQSR.cram
sample=$(basename $bam .sorted.markdup.BQSR.cram)
##提取分裂比对的reads
samtools view -@ 24 -b -F 1294 $bam > $sample.discordants.bam

##提取不正常比对的reads
samtools view -@ 24 -h $bam | \
extractSplitReads_BwaMem -i stdin |samtools sort -@ 24 - > $sample.splitters.bam

##SV鉴定
lumpyexpress \
-B $bam \
-S $sample.splitters.bam \
-D $sample.discordants.bam \
-o $samplel.sv.lumpy.vcf  

##对每个个体SV进行分型
vcftools --vcf $sample.sv.lumpy.vcf \
--indv $sample --recode --recode-INFO-all \
--out $sample  #输出文件前缀
#多个样品可以并行命令，也可以写个循环

svtyper 
-i $sample.recode.vcf \
-B $bam  \
-o $sample.genotype.vcf

#多个样品可以并行命令，也可以写个循环

# 所有样品的vcf合并成群体SV
# 将所有的样品id，列到一个文件里（vcf.list）
# svtools vcfpaste -f vcf.list  > $bam.genotype.vcf
cd ../
echo "##################   end run lumpy   #######################"
