rule lumpy_step1_get_discordants_bam:
    input:
        "{PWD}/data/{id}.bam"
    output:
        "{PWD}/lumpy/{id}.discordants.unsorted.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools view -@ {params.threads} -b -F 1294 {input} > {output}"

rule lumpy_step1_sort_discordants_bam:
    input:
        "{PWD}/lumpy/{id}.discordants.unsorted.bam"
    output:
        "{PWD}/lumpy/{id}.discordants.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools sort -@ {params.threads} {input} -o {output}"

rule lumpy_step2_get_splitters_bam:
    input:
        "{PWD}/data/{id}.bam"
    output:
        "{PWD}/lumpy/{id}.splitters.unsorted.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools view -h {input} | {params.lumpy_path}/extractSplitReads_BwaMem -i stdin | {params.samtools_path}/samtools view -Sb - > {output}"

rule lumpy_step3_sort_splitters_bam:
    input:
        "{PWD}/lumpy/{id}.splitters.unsorted.bam"
    output:
        "{PWD}/lumpy/{id}.splitters.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools sort -@ {params.threads} {input} -o {output}"

rule lumpy_step4_run_lumpy:
    input:
        discordants = "{PWD}/lumpy/{id}.discordants.bam",
        splitters = "{PWD}/lumpy/{id}.splitters.bam",
        raw_bam = "{PWD}/data/{id}.bam"
    output:
        "{PWD}/lumpy/{id}.lumpy.vcf"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"]
    shell:
        """
        {params.lumpy_path}/lumpyexpress -B {input.raw_bam} -S {input.splitters} -D {input.discordants} -o {output}
        """

rule lumpy_genotyped_vcf:
    input:
        vcf = "{PWD}/lumpy/{id}.lumpy.vcf",
        bam = "{PWD}/data/{id}.bam"
    output:
        "{PWD}/lumpy/{id}.genotyped.vcf"
    params:
        env = config["params"]["manta"],
        threads = config["threads"],
        sample = config["samples"]["id"]
    shell:
        """
        vcftools --vcf {input.vcf} \
        --indv {wildcards.PWD}/{params.sample} --recode --recode-INFO-all  \
        --out {wildcards.PWD}/lumpy/{params.sample} && \
        {params.env}/svtyper \
        -i {wildcards.PWD}/lumpy/{params.sample}.recode.vcf \
        -B {input.bam} \
        -o {output}
        """
