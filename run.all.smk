# include the config file
# configfile: "config.yaml"

# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
REF = config['params']['ref_fa']
SAMPLE = config['samples']['id']
PWD = config['params']['PWD']

# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    run_all_files.append("{PWD}/breakdancer/{id}.cfg.SV.output".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/cnvnator/{id}.cnvnator.vcf".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/lumpy/{id}.genotyped.vcf".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/manta/{id}.manta.sv.vcf".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/metasv/{id}.SV.vcf.gz".format(PWD = PWD, id = SAMPLE))

    run_all_files.append("{PWD}/metasv/{id}.manta.gt.vcf".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/metasv/{id}.lumpy.gt.vcf".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/{id}_have_done.txt".format(PWD = PWD, id = SAMPLE))
    run_all_files.append("{PWD}/metasv/{id}.SV.pass.vcf".format(PWD = PWD, id = SAMPLE))
    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input

smk_path = config['params']['smk_path']
include: smk_path+"/breakdancer.smk"
include: smk_path+"/cnvnator.smk"
include: smk_path+"/manta.smk"
include: smk_path+"/metasv.smk"
# include: smk_path+"/pindel.smk"
include: smk_path+"/lumpy.smk"
# include: smk_path+"/annotation.smk"
