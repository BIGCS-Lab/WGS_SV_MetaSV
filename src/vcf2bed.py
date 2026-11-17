import re
import argparse
import gzip
example_text = '''example: python * '''
parser = argparse.ArgumentParser(description=".",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--output_bed', '-o',type= str,help="",required= True,metavar='')
parser.add_argument('--no_complexSV', '-n', action='store_true', default=False, help="no complex SV")
args = parser.parse_args()


def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def vcf2bed(vcf_file, bed_file):
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE =re.compile(r'\bSVTYPE=(<?\w+>?)\b')
    with open(bed_file, 'w') as bedf:
    # with open(vcf_file, 'r') as vcff:
        for line in vcf_file:
            if line.startswith('#'):
                pass
            else:
                chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                svtype = SVTYPE.findall(info)[0].strip("<").strip(">")
                if args.no_complexSV and svtype not in ['DEL', 'DUP', 'INV', 'INS']:
                    continue
                elif '<INS>' in alt or svtype == 'INS':
                    try:
                        svlen = abs(float(SVLEN.findall(info)[0]))
                        end = str(int(pos) + int(svlen))
                    except:
                        end = ENDpos.findall(info)[0]
                else:
                    try:
                        end = ENDpos.findall(info)[0]
                    except:
                        # print(line.strip())
                        svlen = abs(float(SVLEN.findall(info)[0]))
                        end = str(int(pos) + int(svlen))
                print(chrom, int(pos)-1, end, svtype, ".", ".", sep="\t", file=bedf)

if __name__ == '__main__':
    input_vcf = args.input_vcf
    output_bed = args.output_bed
    vcf_file = open_vcf_file(input_vcf)
    vcf2bed(vcf_file, output_bed)
    vcf_file.close()
