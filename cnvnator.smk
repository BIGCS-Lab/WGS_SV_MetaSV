def get_bin_size(depth):
    if 100 >= depth >60:
        bin_size = 30
    elif 50 >= depth > 40:
        bin_size = 50

    elif 40 >= depth > 10:
        bin_size = 100

    elif 10 >= depth > 4:
        bin_size = 500
    else:
        bin_size = 1000
    return bin_size

if str(config["params"]["depth"]).isdigit():
    bin_size = get_bin_size(int(config["params"]["depth"]))
else:
    depth = config["params"]["depth"].strip("x")
    depth = int(depth.strip("X"))
    bin_size = get_bin_size(depth)

# get the alignment information from the bam file
rule cnvnator_step1:
    input:
        bam = "{PWD}/data/{id}.bam",
        ref = REF
    output:
        "{PWD}/{id}/cnvnator/{id}.root"
    params:
        config["params"]["cnvnator"]
    shell:
        "{params}/cnvnator -root {output} -genome {input.ref} -tree {input.bam} -chrom $(seq -f 'chr%g' 1 22) chrX chrY"

rule cnvnator_step2:
    input:
        root = "{PWD}/{id}/cnvnator/{id}.root",
        ref = REF
    output:
        "{PWD}/{id}/cnvnator/{id}.cnvnator"
    params:
        env = config["params"]["cnvnator"],
        bin_size = bin_size
    shell:
        """
        # Step2生成柱形图 
        {params.env}/cnvnator -root {input.root} -his {params.bin_size} -fasta {input.ref} -d {wildcards.PWD}/{wildcards.id}/cnvnator && \
        # Step3统计量计算
        {params.env}/cnvnator -root {input.root} -stat {params.bin_size} -d {wildcards.PWD}/{wildcards.id}/cnvnator  && \
        # Step4 RD信号分割
        {params.env}/cnvnator -root {input.root} -partition {params.bin_size} -ngc -d {wildcards.PWD}/{wildcards.id}/cnvnator  && \
        {params.env}/cnvnator -root {input.root} -call {params.bin_size} -ngc > {output}
        """


rule cnvnator_step3:
    input:
        cnvroot = "{PWD}/{id}/cnvnator/{id}.cnvnator",
        ref = REF
    output:
        "{PWD}/{id}/cnvnator/{id}.cnvnator.vcf"
    params:
        env = config["params"]["cnvnator"],
        # sample = config['samples']['id'],
        chr = config["params"]["ref_fa_chr"]
    shell:
        "{params.env}/perl {params.env}/cnvnator2VCF.pl -prefix {wildcards.PWD}/{wildcards.id} -reference {input.ref} "
        "{input.cnvroot} {params.chr} > {output}"

# rule cnvnator_filter:
#     input:
#         "cnvnator/{id}.cnvnator"
#     output:
#         "cnvnator/{id}.cnvnator.filtered"
#     params:
#         python3 = config["params"]["python3"],
#         smk_path = config["params"]["smk_path"]
#     shell:
#         """
#         {params.python3} {params.smk_path}/src/cnvnator_filter.py {input} {output}
#         """
