
rule manta_step1_configmanta:
    input:
        bam = "{PWD}/data/{id}.bam",
        ref = REF 
    output:
        "{PWD}/manta/{id}.manta.vcf.gz"
    params:
        env = config["params"]['manta'],
        threads = config["threads"]
    shell:
        """
        {params.env}/python {params.env}/configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir {wildcards.PWD}/manta/ && \
        {params.env}/python {wildcards.PWD}/manta/runWorkflow.py -j {params.threads} && cp {wildcards.PWD}/manta/results/variants/diploidSV.vcf.gz {output}
        """

rule manta_convertInversion:
    input:
        vcf = "{PWD}/manta/{id}.manta.vcf.gz",
        ref = REF
    output:
        "{PWD}/manta/{id}.manta.sv.vcf"
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
