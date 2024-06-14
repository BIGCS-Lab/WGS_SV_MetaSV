
rule breakdancer_step1:
    input:
        "{PWD}/data/{id}.bam"
    output:
        "{PWD}/{id}/breakdancer/{id}.cfg"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/perl {params}/bam2cfg.pl {input} > {output}"

rule breakdancer_step2:
    input:
        "{PWD}/{id}/breakdancer/{id}.cfg"
    output:
        "{PWD}/{id}/breakdancer/{id}.cfg.SV.output"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/breakdancer-max {input} > {output}"

