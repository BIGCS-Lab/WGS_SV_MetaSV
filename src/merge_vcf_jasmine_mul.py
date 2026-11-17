#coding=utf-8
from __future__ import print_function 
from __future__ import division
from collections import defaultdict
import argparse
import re
import sys
sys.setrecursionlimit(1000000000) 

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

def line_parser(line):
    end =re.compile(r'\bEND=(-?\d+)\D')
    svlen =re.compile(r'\bSVLEN=(-?\d+)\D')
    lines = line.strip().split('\t')
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *gt_list = lines
    END = end.findall(INFO)[0]
    SVLEN = svlen.findall(INFO)[0]
    return int(END), int(SVLEN), CHROM, int(POS),ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, gt_list


def combine_gt(gt_list, svlen, gt_list2, svlen2):
    if len(gt_list)!= len(gt_list2):
        print("Error: gt_list and gt_list2 have different length")
    c = []
    for i in range(len(gt_list)):
        if gt_list[i] == "0|0:NA:NA:NA:NA" and gt_list2[i] == "0|0:NA:NA:NA:NA":
            c.append("0|0:NA:NA:NA:NA")
        else:
            if gt_list[i] == "0|0:NA:NA:NA:NA":
                c.append(gt_list2[i])
            elif gt_list2[i] == "0|0:NA:NA:NA:NA":
                c.append(gt_list[i])
            elif svlen2 > svlen:
                c.append(gt_list2[i])
            else:
                c.append(gt_list[i])
    # print(c)
    return c

def combine_line(f, baseline, pass_line, cutoff = 0.5):
    try:
        line =  next(f)
    except:
        pass_line.append(baseline.strip())
        return pass_line
    END, SVLEN, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, gt_list = line_parser(baseline)
    END2, SVLEN2, CHROM2,POS2,ID2,REF2,ALT2,QUAL2,FILTER2,INFO2,FORMAT2, gt_list2 = line_parser(line)
    if POS2 >= END or CHROM2 != CHROM or ALT2 != ALT:
        pass_line.append(baseline.strip())
        baseline = line
    elif POS <= POS2 < END:
        if END2 <= END:
            overlap_length = END2 - POS2
            
        elif END <= END2:
            overlap_length = END - POS2
            
        if (cutoff*int(SVLEN2) <= overlap_length) and (cutoff*int(SVLEN) <= overlap_length):
            maxEND = max(END, END2)
            maxSVLEN = max(END, END2) - POS
            new_gt_list = combine_gt(gt_list, SVLEN, gt_list2, SVLEN2)
            new_info = "END=" + str(maxEND) + ";SVLEN=" + str(maxSVLEN) + ";SVTYPE=" + ALT
            new_line = '\t'.join([CHROM, str(POS), ID, REF, ALT, QUAL, FILTER, new_info, FORMAT, '\t'.join(new_gt_list)])
            # pass_line.append(new_line)
            baseline = new_line
        else:
            pass_line.append(baseline.strip())
            baseline = line
    return combine_line(f, baseline, pass_line, cutoff = 0.5)


if __name__ == '__main__':
    args = parser.parse_args()
    input_vcf = args.input_vcf
    overlap_ratio = float(args.overlap_ratio)
    output_vcf = args.output_vcf
    
    with open(input_vcf, 'r') as f:
        pass_line = []
        for line in f:
            if line.startswith('#'):
                pass_line.append(line.strip())
            else:
                baseline = line.strip()
                break
        pass_line = combine_line(f, baseline, pass_line, overlap_ratio)

    with open(output_vcf, 'w') as f:
        for line in pass_line:
            print(line, file=f)