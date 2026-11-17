#coding=utf-8
from __future__ import print_function 
import gzip
import argparse
import re
example_text = '''example:
   python xx.py -i INS.vcf -d DUP.vcf -o INS_DUP.vcf
                '''
parser = argparse.ArgumentParser(description="The script is to merge the INS and DUP vcf into INS DUP vcf files.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--INS_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--DUP_vcf','-d',type=str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')
args = parser.parse_args()

def open_vcf_file(file_name):
    if file_name.endswith('.gz'):
        return gzip.open(file_name, 'rt')
    else:
        return open(file_name, 'r')
    
def cat_vcf_file(INS_vcf,DUP_vcf,output_vcf):
    INS_file = open_vcf_file(INS_vcf)
    DUP_file = open_vcf_file(DUP_vcf)
    output_file = open(output_vcf, 'w')
    for line in INS_file:
        if line.startswith('#'):
            output_file.write(line)
        else:
            output_file.write(line)
    for line in DUP_file:
        if line.startswith('#'):
            pass
        else:
            line = re.sub('DUP', 'INS', line)
            output_file.write(line)
    INS_file.close()
    DUP_file.close()
    output_file.close()

if __name__ == '__main__':
    cat_vcf_file(args.INS_vcf,args.DUP_vcf,args.output_vcf)