
# run annovar software to annotate SVs
rule annnotate_sv_annovar:
    input:
        "metasv_no_pindel/{id}.SV.pass.vcf"
    output:
        "annotation/{id}.SV.pass.annotation.txt"
    params:
        tools = config["params"]["annovar"],
        annovar_db = config["params"]["annovar_db"],
        threads = config["threads"],
        prefix = config["samples"]["id"]
    shell:
        """
        {params.tools}/convert2annovar.pl --format vcf4 {input}  > annotation/{params.prefix}.SV_annovar.input && \
        {params.tools}/annotate_variation.pl --geneanno --neargene 2000  \
        --buildver hg38 --dbtype refGene  --thread {params.threads} \
        --outfile annotation/{params.prefix} \
        --exonsort \
        annotation/{params.prefix}.SV_annovar.input {params.annovar_db} && \
        mv annotation/{params.prefix}.variant_function {output}
        """