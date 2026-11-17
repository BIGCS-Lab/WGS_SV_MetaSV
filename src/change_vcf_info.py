
import sys
import gzip
import re

def args_parse():
    if len(sys.argv)!= 3:
        print("Usage: python change_vcf_info.py input.vcf output.vcf")
        sys.exit(1)

def OpenVcf(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def change_info(inpf, outf):
    ENDpos =re.compile(r'\bEND=(-?\d+)\D')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\D')
    SVTYPE =re.compile(r'\bSVTYPE=(-?\w+)\b')
    for line in inpf:
        if line.startswith('#'):
            outf.write(line)
        else:
            # print(line)
            chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
            try:
                end = ENDpos.findall(info)[0]
            except:
                pass
            try:
                svlen = SVLEN.findall(info)[0]
            except:
                pass
            if ('end' in locals() and'svlen' in locals()):
                pass
            else:
                if 'end' not in locals():
                    end = int(pos) + abs(int(svlen))
                if 'svlen' not in locals():
                    svlen = abs(int(end) - int(pos))
                    
            svtype = SVTYPE.findall(info)[0]
            info = "END=" + str(end) + ";SVLEN=" + str(svlen) + ";SVTYPE=" + svtype
            print(chrom, pos, id, ref, alt, qual, filter, info,*_, sep="\t", file=outf)
            # print(chrom, pos, id, ref, alt, qual, filter, info,*_, sep="\t")
            
if __name__ == '__main__':
    args_parse()
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    inpf = OpenVcf(input_vcf)
    with open(output_vcf, 'w') as outf:
        change_info(inpf, outf)
    inpf.close()
