
rule manta_step1_configmanta:
    input:
        bam = "{PWD}/data/{id}.bam",#.format(config['samples']['id']),
        ref = REF 
    output:
        "{PWD}/{id}/manta/{id}.manta.vcf.gz"
    params:
        config["params"]['manta']
    threads:
        24
    shell:
        """
        {params}/python {params}/configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir {wildcards.PWD}/{wildcards.id}/manta/ && \
        {params}/python {wildcards.PWD}/{wildcards.id}/manta/runWorkflow.py -j 24 && cp {wildcards.PWD}/{wildcards.id}/manta/results/variants/diploidSV.vcf.gz {output}
        """

rule manta_convertInversion:
    input:
        vcf = "{PWD}/{id}/manta/{id}.manta.vcf.gz",
        ref = REF
    output:
        "{PWD}/{id}/manta/{id}.manta.sv.vcf"
    params:
        config["params"]['manta']
    shell:
        """
        {params}/python \
        {params}/convertInversion.py \
        {params}/samtools {input.ref} {input.vcf} > {output}
        """
# rule manta_step2_filter:
#     input:
#         "manta/{id}.manta.vcf.gz"
#     output:
#         "manta/{id}.manta.filtered.vcf"
#     params:
#         python3 = config["params"]['python3'],
#         src = config["params"]['smk_path']
#     shell:
#         """
#         {params.python3} {params.src}/src/manta_filter.py -i {input} -o {output}
#         """
