#coding=utf-8
from __future__ import print_function 
from __future__ import division
import argparse
import os
import glob
import subprocess
from multiprocessing import Pool

example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--path','-p',type=str,help="one or multiple paths that can be separated by commas",required= True,metavar='')
parser.add_argument('--suffix', '-s',type= str,help="",required= True,metavar='')
parser.add_argument('--samtols_tools_path', '-t',action='store_true', default=False, help="samtools tools path")
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


def mul_process(bam_file, read_length):
    input_path = os.path.abspath("input")
    bam_file_path = os.path.abspath(bam_file)
    if not read_length:
        read_length = get_read_length_bam(bam_file)
    else:
        read_length = read_length
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

def make_info_bam(path, suffix, output_file, read_length):
    with open(output_file, 'w') as f:
        bam_files = glob.glob(os.path.join(path, "**", '*.'+suffix), recursive=True)
        os.makedirs("input", exist_ok=True)
        pool = Pool(processes=24)
        all_info = []
        for bam_file in bam_files:
            str = pool.apply_async(mul_process, args=(bam_file, read_length))
            all_info.append(str)
        for info in all_info:
            name, size, read_length = info.get()
            f.write(name + "\t" + size + "\t" + read_length + "\n")
        pool.close()
        pool.join()


if __name__ == '__main__':
    args = parser.parse_args()
    paths = (args.path).strip().split(",")
    suffix = args.suffix
    output_file = args.output_file
    if args.read_length:
        read_length = args.read_length
    else:
        read_length = False
    for path in paths:
        if not os.path.exists(path):
            print("The path {} does not exist".format(path))
        else:
            make_info_bam(path, suffix, output_file, read_length)
