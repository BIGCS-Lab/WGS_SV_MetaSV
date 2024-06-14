import gzip
import argparse
from collections import defaultdict
example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--raw_vcf','-r',type=str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')
parser.add_argument('--metasv_vcf', '-m',type= str,help="",required= True,metavar='')


def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def raw_vcf_info(raw_vcf):
    """
    Parse the raw vcf file and store the information in a dictionary.
    """
    annotations = []
    raw_dict = defaultdict(lambda: defaultdict(list))
    vcff = open_vcf_file(raw_vcf)
    for line in vcff:
        if line.startswith("#"):
            pass
        else:
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample = line.strip().split("\t")
            raw_dict[CHROM][str(POS)].append(FORMAT)
            raw_dict[CHROM][str(POS)].append(sample)
    vcff.close()
    return raw_dict

def replace_GT_info(metasv_vcf, raw_dict, output_vcf):
    """
    replace the GT information in the metasv vcf file with the information from the raw vcf file.
    """
    outf = open(output_vcf, "w")
    vcff = open_vcf_file(metasv_vcf)
    for line in vcff:
        if line.startswith("#"):
            outf.write(line)
        else:
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample = line.strip().split("\t")
            if CHROM in raw_dict:
                if POS in raw_dict[CHROM]:
                    FORMAT2 = raw_dict[CHROM][str(POS)][0].split(":")[0]
                    sample2 = raw_dict[CHROM][str(POS)][1].split(":")[0]
                    print(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT2, sample2, sep="\t", file=outf)
    outf.close()
    vcff.close()
    return "all done"

if __name__ == '__main__':
    args = parser.parse_args()
    raw_dict = raw_vcf_info(args.raw_vcf)
    replace_GT_info(args.metasv_vcf, raw_dict, args.output_vcf)
