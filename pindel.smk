
rule pindel_step1_make_conf:
    input:
        bam = "data/{id}.bam"
    output:
        "pindel/{id}.pindel.conf"
    params:
        sample = config['samples']['id'],
        insert_size = config["params"]["insert_size"]
    shell:
        "echo {input.bam} {params.insert_size} {params.sample} > {output}"

rule pindel_step2_run_pindel:
    input:
        conf = "pindel/{id}.pindel.conf",
        ref = config["params"]["ref_fa"]
    output:
        log = "pindel/{id}.run_pindel.log"
    params:
        env = config["params"]["pindel"],
        sample = config['samples']['id']
    shell:
        "{params.env}/pindel -i {input.conf} -f {input.ref} -o pindel/{params.sample} &> {output.log}"

rule pindel_step3_convert_to_vcf:
    input:
        log = "pindel/{id}.run_pindel.log",
        ref = config["params"]["ref_fa"]
    output:
        "pindel/{id}.pindel.vcf"
    params:
        sample = config['samples']['id'],
        env = config["params"]["pindel"]
    shell:
        "rm {input.log} && {params.env}/pindel2vcf -r {input.ref} -R hg38 -d 20140111 -P pindel/{params.sample} -v {output}"
