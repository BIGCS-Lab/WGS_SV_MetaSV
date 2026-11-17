
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
parser.add_argument('--sv_len', '-l',type= int,help="",metavar='')
# parser.add_argument('--reads_support', '-r',type= int,help="",metavar='')

def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def vcf_filter(vcf_file, out_file, svlen=50):
    PE =re.compile(r'\bPE=(-?\d+)\b') # Number of paired-end reads supporting the variant
    SU =re.compile(r'\bSU=(-?\d+)\b') # Number of pieces of evidence supporting the variant
    SR =re.compile(r'\bSR=(-?\d+)\b') # Number of split reads supporting the variant
    BD =re.compile(r'\bBD=(-?\d+)\b') # Amount of BED evidence supporting the variant
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    with open(out_file, "w") as out_f:
        for line in open_vcf_file(vcf_file):
            if line.startswith("#"):
                out_f.write(line)
            else:
                chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                try:
                    sv_len = SVLEN.findall(info)[0]
                    pe = PE.findall(info)[0]
                    su = SU.findall(info)[0]
                    sr = SR.findall(info)[0]
                    # bd = BD.findall(info)[0]
                except:
                    out_f.write(line)
                    pass
                if abs(int(sv_len)) >= svlen:#  and (int(pe) >= 1 or int(su) >= 1 or int(sr) >= 1):# and int(bd) >= 1:
                    out_f.write(line)

if __name__ == "__main__":
    args = parser.parse_args()
    input_vcf = args.input_vcf
    sv_len = args.sv_len
    # reads_support = args.reads_support
    output_vcf = args.output_vcf
    
    vcf_filter(input_vcf, output_vcf, sv_len)
        