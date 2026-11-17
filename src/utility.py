import re
import gzip

def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def vcf_line_parser(line):
    end =re.compile(r'\bEND=(-?\d+)\D')
    svlen =re.compile(r'\bSVLEN=(-?\d+)\D')
    lines = line.strip().split('\t')
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *gt_list = lines
    END = int(end.findall(INFO)[0])
    SVLEN = abs(int(svlen.findall(INFO)[0]))
    if "INS" in ALT:
        END = int(POS) + int(abs(SVLEN))
    return int(END), int(SVLEN), CHROM, int(POS),ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, gt_list