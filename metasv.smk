from pathlib import Path

def get_cram(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    name = config['samples']['id']
    # cram_files = []
    # for fp in fq_path.iterdir():
    #     # add the fastq file to the files list if it matches
    #     if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".cram"):
    #         print(fp)
    #         cram_files.append(str(fp))
    cram_files = fq_path / f"{name}.sorted.markdup.BQSR.cram"

    if cram_files.exists():
        print(f"file exists: {cram_files}")
    else:
        pattern = f"{name}.*.cram"
        candidates = list(fq_path.glob(pattern))
    
        if len(candidates) == 0:
            raise FileNotFoundError(f"unable to find the cram file of {name} in: {fq_path}")
        elif len(candidates) == 1:
            cram_files = candidates[0]
            print(f"file exists: {cram_files}")
        elif len(candidates) > 1:
            files_str = "\n".join([f"  - {f.name}" for f in candidates])
            raise RuntimeError(f"find more than one cram file of {name}, please check {files_str}")
        else:
            print(f"file not exists: {cram_files}")

    return cram_files 

def get_cram_index(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    name = config['samples']['id']
    # cram_files = []
    # for fp in fq_path.iterdir():
    #     # add the fastq file to the files list if it matches 
    #     if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".crai"):
    #         print(fp)
    #         cram_files.append(str(fp))
    cram_files = fq_path / f"{name}.sorted.markdup.BQSR.cram.crai"

    if cram_files.exists():
        print(f"file exists: {cram_files}")
    else:
        pattern = f"{name}.*.cram.crai"
        candidates = list(fq_path.glob(pattern))
    
        if len(candidates) == 0:
            raise FileNotFoundError(f"unable to find the cram.crai file of {name} in: {fq_path}")
        elif len(candidates) == 1:
            cram_files = candidates[0]
            print(f"file exists: {cram_files}")
        elif len(candidates) > 1:
            files_str = "\n".join([f"  - {f.name}" for f in candidates])
            raise RuntimeError(f"find more than one cram.crai file of {name}, please check {files_str}")
        else:
            print(f"file not exists: {cram_files}")
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
        cram = "{PWD}/data/{id}.cram",
        ref_hg38 = REF
    output:
        "{PWD}/data/{id}.bam"
    threads:
        config['threads']
    params:
        samtools_path = config['params']['samtools']
    shell:
        "{params.samtools_path}/samtools view -@ {threads} -T {input.ref_hg38} -bh {input.cram} -o {output} && {params.samtools_path}/samtools index -@ {threads} {output}"


rule metaSV_merge_vcf:
    input:
        manta_vcf = "{PWD}/manta/{id}.manta.sv.vcf",
        lumpyinput ="{PWD}/lumpy/{id}.lumpy.vcf",
        breakdancerinput= "{PWD}/breakdancer/{id}.cfg.SV.output",
        cnvnatorinput= "{PWD}/cnvnator/{id}.cnvnator.vcf",
        ref = REF
    output:
        metasv_vcf = "{PWD}/metasv/{id}.SV.vcf.gz"
        # breakdancer = "metasv/breakdancer.vcf.gz",
        # manta = "metasv/manta.vcf.gz",
        # cnvnator = "metasv/cnvnator.vcf.gz",
        # lumpy = "metasv/lumpy.vcf.gz"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        sample = config['samples']['id'],
        readlength = config['params']['read_length']
    run:
        shell(
        "{params.env}/python {params.env}/run_metasv.py "
            "--reference {input.ref} --sample {wildcards.PWD}/{params.sample} --disable_assembly --num_threads 2 "
            "--enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} "
            "--outdir {wildcards.PWD}/metasv --workdir {wildcards.PWD}/metasv/tmp_work "
            "--overlap_ratio 0.5 --minsvlen 50 --maxsvlen 10000000 "
            "--breakdancer_native {input.breakdancerinput} "
            "--manta_vcf {input.manta_vcf} "
            "--cnvnator_vcf {input.cnvnatorinput} "
            "--lumpy_vcf {input.lumpyinput} && "
        "mv {wildcards.PWD}/metasv/variants.vcf.gz {output.metasv_vcf} && mv {wildcards.PWD}/metasv/variants.vcf.gz.tbi {output.metasv_vcf}.tbi")

rule metaSV_merge_vcf_filter:
    input:
        "{PWD}/metasv/{id}.SV.vcf.gz"
    output:
        "{PWD}/metasv/{id}.SV.pass.vcf"
    params:
        src = config['params']['smk_path'],
        sample = config['samples']['id'],
        env = config['params']['python3']
    shell:
        """
        {params.env} {params.src}/src/VCF_SeparateFilter.py -i {input} -p {wildcards.PWD}/metasv/{params.sample}_pass -o {output}
        """

rule add_genotype_to_metasv_lumpy:
    input:
        lumpy_vcf = "{PWD}/lumpy/{id}.genotyped.vcf",
        # modify_vcf = "metasv/lumpy.vcf.gz",
        metasv_vcf = "{PWD}/metasv/{id}.SV.vcf.gz"
    output:
        "{PWD}/metasv/{id}.lumpy.gt.vcf"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.lumpy_vcf} -m {wildcards.PWD}/metasv/lumpy.vcf.gz -o {output}"

rule add_genotype_to_metasv_manta:
    input:
        manta_vcf = "{PWD}/manta/{id}.manta.sv.vcf",
        # modify_vcf = "metasv/manta.vcf.gz",
        metasv_vcf = "{PWD}/metasv/{id}.SV.vcf.gz"
    output:
        "{PWD}/metasv/{id}.manta.gt.vcf"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.manta_vcf} -m {wildcards.PWD}/metasv/manta.vcf.gz -o {output}"

rule done:
    input:
        "{PWD}/metasv/{id}.SV.vcf.gz",
        "{PWD}/metasv/{id}.lumpy.gt.vcf",
        "{PWD}/metasv/{id}.manta.gt.vcf"
    output:
        "{PWD}/{id}_have_done.txt"
    params:
        sample = config['samples']['id']
    shell:
        "rm -f {wildcards.PWD}/data/{wildcards.id}.bam* {wildcards.PWD}/cnvnator/{wildcards.id}.root {wildcards.PWD}/lumpy/{wildcards.id}.*.bam && echo {params.sample} have done. > {output}"

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
