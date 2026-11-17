#coding=utf-8
from __future__ import print_function
import argparse
import logging
import os
from collections import defaultdict
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--sample_infomation','-i',type=str,help="the sample information file, colnames: sample_name, depth, read_length,**",required= True,metavar='')
parser.add_argument('--config_yaml', '-c',type= str,help="the config.yaml file template",required= True,metavar='')
parser.add_argument('--run_snakemake_sh', '-r',type= str,help="the run_snakemake.sh file template",required= True,metavar='')
parser.add_argument('--PWD', '-p',type= str,help="the current working directory",required= True,metavar='')
args = parser.parse_args()

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def create_dirs(dirlist):
    if not isinstance(dirlist, list):
        dirlist = [dirlist]
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating directory %s" % (dirname))
            os.makedirs(dirname)
        else:
            logger.info("Directory %s exists" % (dirname))

def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
        logger.info(f"File '{file_path}' has been deleted and will be recreated next time.")
    else:
        logger.info(f"File '{file_path}' does not exist and will be created.")
            
def sample_info_dict(sample_infomation_file):
    """
    convert the sample information file to a dictionary
    """
    samp_info_dict = defaultdict(list)
    with open(sample_infomation_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines = line.strip().split('\t')
                samp_info_dict[lines[0].strip()].extend([lines[1].strip('X'), lines[2].strip('bp')] )
        return samp_info_dict

def change_config_file(path, samplename, depth, read_length, run_snakemake_sh_file, config_yaml_file):
    config_yaml_file_path = os.path.join(path, 'config.yaml')
    _run_name = samplename + '_' + 'run_snakemake.sh'
    run_snakemake_sh_file_path = os.path.join(path, _run_name)
    with open(run_snakemake_sh_file_path, 'w') as outf:
        with open(run_snakemake_sh_file, 'r') as inpf:
            for line in inpf:
                if line.strip(" ").startswith("pwd"):
                    print('pwd='+str(path),file = outf)
                else:
                    outf.write(line)
    with open(config_yaml_file_path, 'w') as outf:
        with open(config_yaml_file, 'r') as configf:
            for line in configf:
                if line.strip(" ").startswith("id"):
                    print('    id: "'+samplename+'"',file = outf)
                elif line.strip(" ").startswith("depth"):
                    print('    depth: '+str(depth),file = outf)
                elif line.strip(" ").startswith("read_length"):
                    print('    read_length: '+str(read_length),file = outf)
                elif line.strip(" ").startswith("PWD"):
                    print('    PWD: '+str(path),file = outf)
                else:
                    outf.write(line)
    

def create_all_dirs_files(samp_info_dict, pwd, config_yaml_file, run_snakemake_sh_file):
    '''
    splitting the script into two parts by depth of bam files, 
    inorder to run the pipeline faster by considering how many jobs took to run in parallel.
    '''
    output_dir = os.path.join(pwd, 'output/01_SV_calling')
    shell_dir = os.path.join(pwd, 'shell')
    dir_7x = os.path.join(shell_dir, '7x')
    dir_30x = os.path.join(shell_dir, '30x')
    create_dirs([dir_7x, dir_30x])
    delete_file(os.path.join(dir_7x, 'all_7x.sh'))
    delete_file(os.path.join(dir_30x, 'all_30x.sh'))
    # create all file in output
    for sample in samp_info_dict:
        sample_dir = os.path.join(output_dir, sample)
        create_dirs([sample_dir])
        depth = samp_info_dict[sample][0].strip('X')
        read_length = samp_info_dict[sample][1].strip('bp')
        change_config_file(sample_dir, sample, depth, read_length, run_snakemake_sh_file, config_yaml_file)
    # create all file in shell
        if 0 < float(depth) <= 10:
            
            with open(os.path.join(dir_7x, 'all_7x.sh'), 'a') as outf:
                _run_name = os.path.join(sample_dir, str(sample) + '_' + 'run_snakemake.sh')
                print(f'sh {_run_name}', file = outf)
        elif 10 < float(depth):
            
            with open(os.path.join(dir_30x, 'all_30x.sh'), 'a') as outf:
                _run_name = os.path.join(sample_dir, str(sample) + '_' + 'run_snakemake.sh')
                print(f'sh {_run_name}', file = outf)
            
if __name__ == '__main__':
    sample_infomation_file = args.sample_infomation
    config_yaml_file = args.config_yaml
    run_snakemake_sh_file = args.run_snakemake_sh
    pwd = os.path.abspath(args.PWD)
    samp_info_dict = sample_info_dict(sample_infomation_file)
    create_all_dirs_files(samp_info_dict, pwd, config_yaml_file, run_snakemake_sh_file)