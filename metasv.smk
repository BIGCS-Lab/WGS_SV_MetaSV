from pathlib import Path

def get_cram(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    # for name in config['samples']['id']:
    # name = config['samples']['id']
    name = wildcards.id
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".cram"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files 

def get_cram_index(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    # for name in config['samples']['id']:
    # name = config['samples']['id']
    name = wildcards.id
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches 
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".crai"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files

# link the cram and cram index files to the data directory
rule link_cram:
    input:
        get_cram
    output:
        "{PWD}/data/{id}.cram"
    run:
        shell("ln -s {input} {output}")

rule link_cram_index:
    input:
        get_cram_index
    output:
        "{PWD}/data/{id}.cram.crai"
    run:
        shell("ln -s {input} {output}")

# convert the cram file to bam file
rule bam_from_cram:
    input:
        "{PWD}/data/{id}.cram"
    output:
        "{PWD}/data/{id}.bam"
    threads:
        config['threads']
    params:
        samtools_path = config['params']['samtools']
    shell:
        "{params.samtools_path}/samtools view -@ {threads} -bh {input} -o {output} && {params.samtools_path}/samtools index -@ {threads} {output}"


rule metaSV_merge_vcf:
    input:
        manta_vcf = "{PWD}/{id}/manta/{id}.manta.sv.vcf",
        lumpyinput ="{PWD}/{id}/lumpy/{id}.lumpy.vcf",
        breakdancerinput= "{PWD}/{id}/breakdancer/{id}.cfg.SV.output",
        cnvnatorinput= "{PWD}/{id}/cnvnator/{id}.cnvnator.vcf",
        ref = REF
    output:
        metasv_vcf = "{PWD}/{id}/metasv/{id}.SV.vcf.gz"
        # breakdancer = "metasv/breakdancer.vcf.gz",
        # manta = "metasv/manta.vcf.gz",
        # cnvnator = "metasv/cnvnator.vcf.gz",
        # lumpy = "metasv/lumpy.vcf.gz"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        # sample = config['samples']['id'],
        readlength = config['params']['read_length']
    run:
        shell(
        "{params.env}/python {params.env}/run_metasv.py "
            "--reference {input.ref} --sample {wildcards.PWD}/{wildcards.id} --disable_assembly --num_threads {params.threads} "
            "--enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} "
            "--outdir {wildcards.PWD}/{wildcards.id}/metasv --workdir {wildcards.PWD}/{wildcards.id}/metasv/tmp_work "
            "--overlap_ratio 0.5 --minsvlen 50 --maxsvlen 10000000 "
            "--breakdancer_native {input.breakdancerinput} "
            "--manta_vcf {input.manta_vcf} "
            "--cnvnator_vcf {input.cnvnatorinput} "
            "--lumpy_vcf {input.lumpyinput} && "
        "mv {wildcards.PWD}/{wildcards.id}/metasv/variants.vcf.gz {output.metasv_vcf} && mv {wildcards.PWD}/{wildcards.id}/metasv/variants.vcf.gz.tbi {output.metasv_vcf}.tbi")

rule add_genotype_to_metasv_lumpy:
    input:
        lumpy_vcf = "{PWD}/{id}/lumpy/{id}.genotyped.vcf",
        # modify_vcf = "metasv/lumpy.vcf.gz",
        metasv_vcf = "{PWD}/{id}/metasv/{id}.SV.vcf.gz"
    output:
        "{PWD}/{id}/metasv/{id}.lumpy.gt.vcf"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.lumpy_vcf} -m {wildcards.PWD}/{wildcards.id}/metasv/lumpy.vcf.gz -o {output}"

rule add_genotype_to_metasv_manta:
    input:
        manta_vcf = "{PWD}/{id}/manta/{id}.manta.sv.vcf",
        # modify_vcf = "metasv/manta.vcf.gz",
        metasv_vcf = "{PWD}/{id}/metasv/{id}.SV.vcf.gz"
    output:
        "{PWD}/{id}/metasv/{id}.manta.gt.vcf"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.manta_vcf} -m {wildcards.PWD}/{wildcards.id}/metasv/manta.vcf.gz -o {output}"


rule rm_bam_done:
    input:
        "{PWD}/{id}/metasv/{id}.SV.vcf.gz",
        "{PWD}/{id}/metasv/{id}.lumpy.gt.vcf",
        "{PWD}/{id}/metasv/{id}.manta.gt.vcf"
    output:
        "{PWD}/{id}/{id}_have_done.txt"
    # params:
        # sample = config['samples']['id']
    shell:
        "rm -f {wildcards.PWD}/{wildcards.id}/lumpy/{wildcards.id}.*.bam* {wildcards.PWD}/data/{wildcards.id}.bam* {wildcards.PWD}/{wildcards.id}/cnvnator/{wildcards.id}.root && echo {wildcards.id} have done. > {output}"


# rule split_vcf_by_svtype:
#     input:
#         "metasv/{id}.metasv.genotype.vcf".format(id = config['samples']['id'])
#     output:
#         expand("metasv/{id}_{svtype}.vcf", id = config['samples']['id'], svtype=["DEL", "DUP", "INV", "INS"])
#     params:
#         sample = config['samples']['id'],
#         python = config['params']['python3'],
#         src = config['params']['smk_path']
#     run:
#         shell(
#             "{params.python} {params.src}/src/split_vcf_bysvtype.py {input} metasv/{params.sample}")
