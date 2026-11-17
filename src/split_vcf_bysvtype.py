import sys
import gzip
from collections import defaultdict
import re

def open_vcf_file(filename):
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def svtype_dict(vcf_f):
    # SVTYPE = re.compile(r'(?<=\bSVTYPE=)(?P<svtype>[\w:\-_]+)(?=;)')
    SVTYPE = re.compile(r'\bSVTYPE=(-?\w+)\b')
    SV_dict = defaultdict(list)
    vcf_header = []
    for line in vcf_f:
        if line.startswith("#"):
            vcf_header.append(line)
        else:
            chr, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split("\t")
            svtype = SVTYPE.findall(info)[0]
            SV_dict[svtype].append(line)
    return SV_dict, vcf_header

def write_svtype_vcf(SV_dict, vcf_header, outprefix):
    import os
    if os.path.exists(outprefix+"_stat.txt"):
        os.remove(outprefix+"_stat.txt")
    for svtype in SV_dict:
        with open(outprefix+"_"+svtype+".vcf", "w") as out_f:
            out_f.writelines(vcf_header)
            out_f.writelines(SV_dict[svtype])
        with open(outprefix+"_stat.txt", "a") as out_stat:
            out_stat.write(svtype + "\t" + str(len(SV_dict[svtype])) + "\n")
        
def main():
    vcf_file = sys.argv[1]
    outprefix = sys.argv[2]
    vcf_f = open_vcf_file(vcf_file)
    SV_dict, vcf_header = svtype_dict(vcf_f)
    write_svtype_vcf(SV_dict, vcf_header, outprefix)
    vcf_f.close()

if __name__ == "__main__":
    main()
