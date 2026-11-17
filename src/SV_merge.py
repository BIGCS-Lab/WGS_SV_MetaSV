
import gzip
from collections import defaultdict
import re


def open_vcf_file(filename):
    _vcf_f = None

    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "r")
    else:
        _vcf_f = open(filename, "r")

    return _vcf_f

def read_header(vcf_f):
    header = []
    for line in vcf_f:
        if line.startswith("#"):
            header.append(line)
    return header

def multi_vcf_dict(vcf_string):
    # get header
    for vcf_file in vcf_string.split(",")[0]:
        vcf_f = open_vcf_file(vcf_file)
        header = read_header(vcf_f)
    vcf_f.close()
    
    # convert multi-vcf to a big dictionary
    vcf_dict = defaultdict(dict)
    ENDpos =re.compile(r'\bEND=(-?\d+)\D')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\D')
    for vcf_file in vcf_string.split(","):
        vcf_f = open_vcf_file(vcf_file)
        for line in vcf_f:
            if not line.startswith("#"):
                sv_len = SVLEN.search(info)
                chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                vcf_dict[chrom][int(pos)] = line.strip()

