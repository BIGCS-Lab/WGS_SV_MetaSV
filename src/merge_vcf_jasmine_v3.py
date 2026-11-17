#coding=utf-8
from __future__ import print_function 
from __future__ import division
from collections import defaultdict
import argparse
import re
import sys
# sys.setrecursionlimit(1000000000) 

example_text = '''example:
   python *
                '''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--overlap_ratio', '-r',type= str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')

def open_vcf_file(filename):
    import gzip
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def line_parser(line):
    end =re.compile(r'\bEND=(-?\d+)\D')
    svlen =re.compile(r'\bSVLEN=(-?\d+)\D')
    svtype = re.compile(r'\bSVTYPE=(-?\w+)\b')
    lines = line.strip().split('\t')
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *gt_list = lines
    END = int(end.findall(INFO)[0])
    SVLEN = abs(int(svlen.findall(INFO)[0]))
    SVTYPE = svtype.findall(INFO)[0]
    if ("INS" or "BND" or "TRA") in ALT:
        if not END:
            END = int(POS) + int(abs(SVLEN))
    return int(END), int(SVLEN), SVTYPE, CHROM, int(POS),ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, gt_list


def combine_gt(gt_list, svlen, gt_list2, svlen2):
    if len(gt_list)!= len(gt_list2):
        print("Error: gt_list and gt_list2 have different length")
    c = []
    for i in range(len(gt_list)):
        if gt_list[i] == gt_list2[i]:
            c.append(gt_list[i])
        else:
            if gt_list[i].split(':')[0] == "0/0":
                c.append(gt_list2[i])
            elif gt_list2[i].split(':')[0] == "0/0":
                c.append(gt_list[i])
            elif svlen2 > svlen:
                c.append(gt_list2[i])
            else:
                c.append(gt_list[i])
    # print(c)
    return c

def two_loci_overlap(POS, END, POS2, END2, cutoff):
    if POS < END2 and END > POS2:
        overlap_length = min(END, END2) - max(POS, POS2)
        loci_1_length = END - POS
        loci_2_length = END2 - POS2
        if loci_1_length == 0 or loci_2_length == 0:
            return False
        elif overlap_length / max(loci_1_length, loci_2_length) >= float(cutoff):
            return True
    else:
        return False


def combine_line(baseline, line, cutoff, outf):
    # try:
    #     line =  next(f)
    # except:
    #     # pass_line.append(baseline.strip())
    #     # return pass_line
    #     print(baseline.strip(), file=outf)
    END, SVLEN, SVTYPE, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, gt_list = line_parser(baseline)
    END2, SVLEN2, SVTYPE2, CHROM2,POS2,ID2,REF2,ALT2,QUAL2,FILTER2,INFO2,FORMAT2, gt_list2 = line_parser(line)
    
    if (CHROM2 != CHROM) or (SVTYPE2 != SVTYPE) or (("INS" or "BND" or "TRA") in ALT) or (("INS" or "BND" or "TRA") in SVTYPE):
        INS_info = "END=" + str(END) + ";SVLEN=" + str(SVLEN) + ";SVTYPE=" + SVTYPE
        baseline = '\t'.join([CHROM, str(POS), ID, REF, ALT, QUAL, FILTER, INS_info, FORMAT, '\t'.join(gt_list)])
        # pass_line.append(baseline.strip())
        print(baseline.strip(), file=outf)
        baseline = line
        return baseline
    elif two_loci_overlap(POS, END, POS2, END2, cutoff):
        maxEND = max(END, END2)
        maxSVLEN = max(END, END2) - POS
        new_gt_list = combine_gt(gt_list, SVLEN, gt_list2, SVLEN2)
        new_info = "END=" + str(maxEND) + ";SVLEN=" + str(maxSVLEN) + ";SVTYPE=" + SVTYPE
        new_line = '\t'.join([CHROM, str(POS), ID, REF, ALT, QUAL, FILTER, new_info, FORMAT, '\t'.join(new_gt_list)])
        # pass_line.append(new_line)
        baseline = new_line
        return baseline
    else:
        INS_info = "END=" + str(END) + ";SVLEN=" + str(SVLEN) + ";SVTYPE=" + SVTYPE
        baseline = '\t'.join([CHROM, str(POS), ID, REF, ALT, QUAL, FILTER, INS_info, FORMAT, '\t'.join(gt_list)])
        # pass_line.append(baseline.strip())
        print(baseline.strip(), file=outf)
        baseline = line
        return baseline
    # return combine_line(f, baseline, pass_line, cutoff, outf)


if __name__ == '__main__':
    args = parser.parse_args()
    input_vcf = args.input_vcf
    overlap_ratio = float(args.overlap_ratio)
    output_vcf = args.output_vcf
    
    f = open_vcf_file(input_vcf)
    with open(output_vcf, 'w') as outf:
        for i, line in enumerate(f):
            if i == 0:
                baseline = line
            else:
                baseline = combine_line(baseline, line, overlap_ratio, outf)
        # pass_line = combine_line(f, baseline, pass_line, overlap_ratio)
        # combine_line(f, baseline, pass_line, overlap_ratio, outf)
    f.close()

    # with open(output_vcf, 'w') as f:
    #     for line in pass_line:
    #         print(line, file=f)