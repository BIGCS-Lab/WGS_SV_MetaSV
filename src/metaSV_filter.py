
import gzip
import re
import sys
import argparse
example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')
# parser.add_argument('--reads_support', '-r',type= int,help="",metavar='')

def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def vcf_filter(vcf_file, out_file):
    with open(out_file, "w") as out_f:
        for line in open_vcf_file(vcf_file):
            if line.startswith("#"):
                out_f.write(line)
            else:
                if ("BreakDancer" in line) and ("Lumpy" in line) and ("CNVnator" not in line):
                    pass
                else:
                    out_f.write(line)

if __name__ == "__main__":
    args = parser.parse_args()
    input_vcf = args.input_vcf
    # sv_len = args.sv_len
    # reads_support = args.reads_support
    output_vcf = args.output_vcf
    
    vcf_filter(input_vcf, output_vcf)