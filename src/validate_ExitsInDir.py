import os
import shutil

import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--version', action='version', version='%(prog)s 20240702')
parser.add_argument('--path','-p',type=str,help="the path of the log files",required= True,metavar='')
parser.add_argument('--all_sample_id_file','-a',type=str,help="the file of all sample id",required= True,metavar='')
parser.add_argument('--prefix','-p',type=str,help="the prefix of the output files",required= True,metavar='')
args = parser.parse_args()


def validate_done_file(directory_path, file_name):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):
        return True
    else:
        return False
    
def delete_directory(directory_path):
    """
    delete directory recursively
    """
    try:
        shutil.rmtree(directory_path)
        print(f'warning: you have rm the {directory_path} directory.')
    except Exception as e:
        pass
    
def rm_undone_Data_Manta_dir(path, sampleid):
    done_list = []
    undone_list = []
    breakdancer_dir = os.path.join(path, "breakdancer")
    cnvnator_dir = os.path.join(path, "cnvnator")
    lumpy_dir = os.path.join(path, "lumpy")
    manta_dir = os.path.join(path, "manta")
    metasv_dir = os.path.join(path, "metasv")

    breakdancer_file = sampleid + ".cfg.SV.output"
    cnvnator_file = sampleid + ".cnvnator.vcf"
    lumpy_file = sampleid + ".genotyped.vcf"
    manta_file = sampleid + ".manta.vcf.gz"
    metasv_file = sampleid + ".SV.pass.vcf"
    
    if (validate_done_file(breakdancer_dir, breakdancer_file) and 
        validate_done_file(cnvnator_dir, cnvnator_file) and 
        validate_done_file(lumpy_dir, lumpy_file) and 
        validate_done_file(metasv_dir, metasv_file) and
        validate_done_file(manta_dir, manta_file)):
        done_list.append(sampleid)
    else:
        # print(f"{sampleid} is undone in {path}")
        undone_list.append(sampleid)
    if validate_done_file(manta_dir, manta_file):
        pass
    else:
        delete_directory(os.path.join(path, 'manta'))
    return done_list, undone_list

if __name__ == '__main__':
    path = args.path
    all_sample_id_file = args.all_sample_id_file
    with open(all_sample_id_file, 'r') as f:
        for line in f:
            sampleid = line.strip().split('\t')[0]
            path_sampleid = os.path.join(path, sampleid)
            done_list, undone_list = rm_undone_Data_Manta_dir(path_sampleid, sampleid)
            with open(f'{args.prefix}_done_list.txt', 'a') as f1:
                print(*done_list, sep='\n', file=f1)
            with open(f'{args.prefix}_undone_list.txt', 'a') as f2:
                print(*undone_list, sep='\n', file=f2)


















