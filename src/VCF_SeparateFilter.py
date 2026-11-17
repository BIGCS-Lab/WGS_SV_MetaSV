#coding=utf-8
from __future__ import print_function 
from __future__ import division
from collections import defaultdict
import gzip
import argparse
import re
example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is to filter SVs with PASS and separete SV types from VCF file.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')
args = parser.parse_args()

input_vcf = args.input_vcf
prefix = args.prefix
output_vcf = args.output_vcf

# prefix = output.split(".")[0]
SV_len =re.compile(r'\bSVLEN=(-?\d+)')
annotation = []
SVtypes_dict = defaultdict(list)
SVtypes_stat = defaultdict(dict)
sv_len_dict = defaultdict(list)

CNV_outf = open(prefix + '.CNV.vcf', 'w')
with open(output_vcf, 'w') as outf:
    with gzip.open(input_vcf, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if line.startswith('#'):
                annotation.append(line)
                outf.write(line)
                CNV_outf.write(line)
            else:
                sv_len = abs(int(SV_len.findall(line.strip())[0]))
                if "PASS" in fields[6]:
                    outf.write(line)
                    SVtypes_dict[fields[4]].append(line)
                    SVtypes_stat[fields[4]] = SVtypes_stat.get(fields[4], {'count': 0, 'len': 0})
                    SVtypes_stat[fields[4]]['count'] += 1
                    SVtypes_stat[fields[4]]['len'] += abs(int(sv_len))
                    sv_len_dict[fields[4]].append(abs(int(sv_len)))
                    if "CNVnator" in line:
                        CNV_outf.write(line)
                    
                # retain some results for specific SV types that're not PASS
                elif (10000000 > sv_len > 50) and (("<INS>" in line and "Manta" in line) or ("<INV>" in line and "Manta" in line)):
                    """
                    using Mantas INS and INV results to take advantage of its strength in detecting INS and INV
                    """
                    outf.write(line)
                    SVtypes_dict[fields[4]].append(line)
                    SVtypes_stat[fields[4]] = SVtypes_stat.get(fields[4], {'count': 0, 'len': 0})
                    SVtypes_stat[fields[4]]['count'] += 1
                    SVtypes_stat[fields[4]]['len'] += abs(int(sv_len))
                    sv_len_dict[fields[4]].append(abs(int(sv_len)))
                elif (10000000 > sv_len > 10000) and "CNVnator" in line:
                    """
                    retaining the results for CNVnator greater than 10 kb and lesser than 10 Mb
                    """
                    outf.write(line)
                    CNV_outf.write(line)
                    SVtypes_dict[fields[4]].append(line)
                    SVtypes_stat[fields[4]] = SVtypes_stat.get(fields[4], {'count': 0, 'len': 0})
                    SVtypes_stat[fields[4]]['count'] += 1
                    SVtypes_stat[fields[4]]['len'] += abs(int(sv_len))
                    sv_len_dict[fields[4]].append(abs(int(sv_len)))
                
CNV_outf.close()
# convert CNV vcf file to bed file
def vcf2bed(vcf_file, bed_file):
    ENDpos =re.compile(r'\bEND=(-?\d+)\D')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\D')
    with open(bed_file, 'w') as bedf:
        with open(vcf_file, 'r') as vcff:
            for line in vcff:
                if line.startswith('#'):
                    pass
                else:
                    chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                    if '<INS>' in alt:
                        sv_len = SVLEN.findall(info)[0]
                        end = str(int(pos) + int(sv_len))
                    else:
                        end = ENDpos.findall(info)[0]
                    print(chrom, int(pos)-1, end, alt, ".", ".", sep="\t", file=bedf)
vcf2bed(prefix + '.CNV.vcf', prefix + '.CNV.bed')

for SVtype in SVtypes_dict:
    with open(prefix + "." + SVtype.strip("<").strip(">") + ".vcf", 'w') as outf:
        outf.writelines(annotation)
        outf.writelines(SVtypes_dict[SVtype])
    with open(prefix + "." + SVtype.strip("<").strip(">") + ".bed", 'w') as bedf:
        ENDpos =re.compile(r'(?<=\bEND=)(?P<end>\d+)')
        SVLEN =re.compile(r'\bSVLEN=(-?\d+)\D')
        for line in SVtypes_dict[SVtype]:
            if line.startswith('#'):
                pass
            else:
                chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                try:
                    end = ENDpos.findall(info)[0]
                except:
                    sv_len = SVLEN.findall(info)[0]
                    end = str(int(pos) + abs(int(sv_len)))
                print(chrom, int(pos)-1, end, alt, ".", ".", sep="\t", file=bedf)

with open(prefix + '_svtype_stats.txt', 'w') as f:
    print("SVtype", "Count", "Length", sep="\t", file=f)
    for svtype in SVtypes_stat:
        print(svtype.strip(">").strip("<"), SVtypes_stat[svtype]['count'], SVtypes_stat[svtype]['len'], file=f, sep="\t")
        
with open(prefix + '_svtype_len.txt', 'w') as f:
    for i in sv_len_dict:
        for j in sv_len_dict[i]:
            print(i.strip(">").strip("<"), j, prefix, file=f, sep="\t")
            
