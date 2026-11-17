#coding=utf-8
from __future__ import print_function 
from __future__ import division
import argparse
import os
import glob
import subprocess

example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--path','-p',type=str,help="",required= True,metavar='')
parser.add_argument('--suffix', '-s',type= str,help="",required= True,metavar='')
parser.add_argument('--samtols_tools_paht', '-t',action='store_true', default=False, help="samtools tools path")
parser.add_argument('--read_length', '-l', type=int, help="the read length of the bam file")
parser.add_argument('--output_file', '-o',type= str,help="",required= True,metavar='')

def get_read_length_bam(bam_file, samtools = "samtools"):
    cmd = samtools + " view " + bam_file + " |less|head -n 100"
    _res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    # import re
    lens = []
    for line in _res.stdout.strip().split('\n'):
        read_length = len(line.strip().split('\t')[9])
        lens.append(read_length)
    return max(lens)

def make_info_bam(path, suffix, output_file, read_length):
    with open(output_file, 'w') as f:
        bam_files = glob.glob(os.path.join(path, "**", '*.'+suffix), recursive=True)
        for bam_file in bam_files:
            os.makedirs("input", exist_ok=True)
            input_path = os.path.abspath("input")
            bam_file_path = os.path.abspath(bam_file)
            if not read_length:
                read_length2 = get_read_length_bam(bam_file)
            else:
                read_length2 = read_length
            if not os.path.exists(os.path.join(input_path, os.path.basename(bam_file))):
                os.symlink(bam_file_path, os.path.join(input_path, os.path.basename(bam_file)))
            if not os.path.exists(os.path.join(input_path, os.path.basename(bam_file)+".crai")):
                os.symlink(bam_file_path+".crai", os.path.join(input_path, os.path.basename(bam_file)+".crai"))
            file_size_bytes = os.path.getsize(bam_file)
            size = file_size_bytes / (1024 * 1024 * 1024)
            name = os.path.basename(bam_file).split(".")[0]
            if 0 < float(size) < 10:
                return name, "7X", str(read_length)
            elif 10 < float(size) < 20:
                return name, "15X", str(read_length)
            elif 20 < float(size):
                return name, "30X", str(read_length)


if __name__ == '__main__':
    args = parser.parse_args()
    path = args.path
    suffix = args.suffix
    output_file = args.output_file
    if args.read_length:
        read_length = args.read_length
    else:
        read_length = False
    make_info_bam(path, suffix, output_file, read_length)
