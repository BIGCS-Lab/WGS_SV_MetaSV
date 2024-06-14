#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 
############## Step1  从bam文件中提取比对信息
bam=../00116021276M31BFF2A.sorted.markdup.BQSR.cram
sample=$(basename $bam .sorted.markdup.BQSR.cram)
ref=/WORK/gzfezx_shhli_3/BioDatahub/human_reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
mkdir cnvnator
cd cnvnator

cnvnator=/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/cnvnator
tools_dir=/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin
$cnvnator -root $sample.root -tree $bam  -chrom $(seq -f 'chr%g' 1 22) chrX chrY

# -root out.root-指定输出ROOT文件。
# -chrom name1 ...-指定染色体名称。
# -tree file1.bam ...-指定bam文件的名称。
echo '------------------------------------------------------'
echo 
echo "############## Step2生成柱形图 "
# cnvnator -root file.root -his 500 -d dir_with_genome_fa/  
$cnvnator -root $sample.root -his 100 -fasta $ref
# cnvnator -root file.root -his 500 -chrom 1 2 3 4 -fasta file_genome.fa.gz
# -his 参数的值，也就是bin_size的值，这个可以根据CNVnator作者的建议
# 根据测序深度的不同进行调整：100X: 30；20-30X: 100；4-6X: 500；再低就1000。
# 如果使用 -d 那就后面接参考基因组所在目录。
# 如果使用 -fasta 参数就在后面跟参考基因组的文件
echo '------------------------------------------------------'
echo 
echo "############# Step3统计量计算"
$cnvnator -root $sample.root -stat 100 
# 必须先完成此步骤，然后再进行分区和CNV调用。
echo '------------------------------------------------------'
echo 
echo "############# Step4 RD信号分割"
$cnvnator -root $sample.root -partition 100 -ngc
# 选项-ngc指定不使用GC校正的RD信号。分区是最耗时的步骤。
echo '------------------------------------------------------'
echo 
echo "############# Step5 call CNV"
$cnvnator -root $sample.root -call 100 -ngc > $sample.cnvnator
# 结果会输出到sample.cnvnator里。
# CNV_type， coordinates ，CNV_size， normalized_RD， e-val1， e-val2， e-val3 ，e-val4 ，q0.
echo '------------------------------------------------------'
echo 
echo "还可以把sample.cnvnator文件 转换成vcf文件"
$tools_dir/cnvnator2VCF.pl -prefix $sample -reference $ref $sample.cnvnator > $sample.vcf
# -prefix指定要附加到输出VCF中ID字段的前缀字符。
# reference代表您使用的参考基因组的名称，例如GRCh37，hg19等。
# sample.cnvnator是带有CNV调用的CNVnator输出文件
