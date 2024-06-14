# include the config file
# configfile: "config.yaml"

# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
REF = config['params']['ref_fa']
SAMPLE_list = config['samples']['id_list']
PWD = config['params']['PWD']

# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    run_all_files.extend([f'{PWD}/{id}/metasv/{id}.manta.gt.vcf' for id in SAMPLE_list])
    run_all_files.extend([f'{PWD}/{id}/metasv/{id}.lumpy.gt.vcf' for id in SAMPLE_list])
    run_all_files.extend([f'{PWD}/{id}/metasv/{id}.SV.vcf.gz' for id in SAMPLE_list])
    run_all_files.extend([f'{PWD}/{id}/{id}_have_done.txt' for id in SAMPLE_list])
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
