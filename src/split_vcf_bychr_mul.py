#coding=utf-8
from __future__ import print_function 
from __future__ import division
import os
import argparse
from multiprocessing import Pool

def parse_args():
    parser = argparse.ArgumentParser(description='Split VCF by chromosome')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    return parser.parse_args()

def mul_process(line, chr_list):
    if line.startswith('#'):
        for chr in chr_list:
            filef = os.path.join('chrs', chr + '.vcf')
            with open(filef, 'a') as f:
                f.write(line)
    line_list = line.strip().split('\t')
    chr = line_list[0]
    if chr in chr_list:
        filef = os.path.join('chrs', chr + '.vcf')
        with open(filef, 'a') as f:
            f.write(line)

def split_vcf_bychr(input_file, chr_list):
    with open(input_file, 'r') as f:
        pool = Pool(processes=24)
        for line in f:
            pool.apply_async(mul_process, args=(line, chr_list), callback= setcallback)
        pool.close()
        pool.join()

if __name__ == '__main__':
    args = parse_args()
    chr_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    os.makedirs('chrs', exist_ok=True)
    split_vcf_bychr(args.input, chr_list)

