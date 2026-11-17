#coding=utf-8
from __future__ import print_function
import argparse
import gzip
import os
from multiprocessing import Pool

# from collections import defaultdict
import re
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--vcffile','-v',type=str,help="the sv vcf file",required= True,metavar='')
parser.add_argument('--outputpath','-p',type=str,help="the path of output", metavar='')
parser.add_argument('--onlychangeInfo','-c',type=str,help="only change the info of sv", metavar='')

args = parser.parse_args()

def OpenVcf(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def modify_gtlist(gtlist):
    for i in range(len(gtlist)):
        gtlist[i] = gtlist[i].strip().split(':')[0]
        if gtlist[i] == './.':
            gtlist[i] = '0/0'
        elif gtlist[i] == '.|.':
            gtlist[i] = '0/0'
    return gtlist


def change_info(line, header, path, chr_list):
    # chr_list = [''.join(chr(i) for i in range(ord('A'), ord('Z')+1)), 'X', 'Y']
    # chr_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    # chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE =re.compile(r'\bSVTYPE=(-?\w+)\b')
    # for line in inpf:
    if line.startswith('#'):
        with open(os.path.join(path, header), 'a') as f:
            f.write(line)
    else:
        # print(line)
        chrom, pos, id, ref, alt, qual, filter, info, *gtlist = line.strip().split('\t')
        if chrom in chr_list:
            try:
                end = ENDpos.findall(info)[0]
            except:
                # continue
                pass
                # end = "None"
            try:
                svlen = SVLEN.findall(info)[0]
                end = ENDpos.findall(info)[0]
            except:
                svlen = int(end) - int(pos)
            svtype = SVTYPE.findall(info)[0]
            gt = modify_gtlist(gtlist)
            info = "END=" + end + ";SVLEN=" + str(svlen) + ";SVTYPE=" + svtype
            # print(chrom, pos, id, ref, alt, qual, filter, info,*_, sep="\t", file=outf)
            # return svtype, chrom, pos, id, ref, alt, qual, filter, info,gt
            _fold = os.path.join(path, chrom)
            chr_sv_file = os.path.join(_fold, chrom + '_' + svtype + '.vcf')
            with open(chr_sv_file, 'a') as f:
                print(chrom, pos, id, ref, alt, qual, filter, info,*gt, sep="\t", file=f)
                
def only_change_info(line, outf, headerf, chr_list):
    # chr_list = [''.join(chr(i) for i in range(ord('A'), ord('Z')+1)), 'X', 'Y']
    # chr_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    # chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE =re.compile(r'\bSVTYPE=(-?\w+)\b')
    # for line in inpf:
    if line.startswith('##'):
        print(line.strip(), file=headerf)
    elif line.startswith('#CHROM'):
        chrom, pos, id, ref, alt, qual, filter, info, *gtlist = line.strip().split('\t')
        gt = [os.path.basename(x) for x in gtlist]
        print(chrom, pos, id, ref, alt, qual, filter, info, *gt, sep="\t", file=headerf)
    else:
        # print(line)
        chrom, pos, id, ref, alt, qual, filter, info, *gtlist = line.strip().split('\t')
        if chrom in chr_list:
            try:
                end = ENDpos.findall(info)[0]
            except:
                # continue
                pass
                # end = "None"
            try:
                svlen = SVLEN.findall(info)[0]
                end = ENDpos.findall(info)[0]
            except:
                svlen = int(end) - int(pos)
            svtype = SVTYPE.findall(info)[0]
            gt = modify_gtlist(gtlist)
            info = "END=" + end + ";SVLEN=" + str(svlen) + ";SVTYPE=" + svtype
            # print(chrom, pos, id, ref, alt, qual, filter, info,*_, sep="\t", file=outf)
            # return svtype, chrom, pos, id, ref, alt, qual, filter, info,gt
            # _fold = os.path.join(path, chrom)
            # chr_sv_file = os.path.join(_fold, chrom + '_' + svtype + '.vcf')
            # with open(chr_sv_file, 'a') as f:
            print(chrom, pos, id, ref, alt, qual, filter, info,*gt, sep="\t", file=outf)

def only_change(vcf_file, output, chr_list=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']):
    vcffile = OpenVcf(vcf_file)
    # pool = Pool(processes=24)
    # for line in vcffile:
    #     pool.apply_async(only_change_info, args=(line, chr_list))
    # pool.close()
    # pool.join()
    with open(output, 'w') as outf:
        with open(output+'.header', 'w') as headerf:
            for line in vcffile:
                only_change_info(line, outf, headerf, chr_list)
    vcffile.close()

def sv_class(vcf_file, header='hearder.txt', path='./', chr_list=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']):
    vcffile = OpenVcf(vcf_file)
    if not os.path.exists(path):
        os.mkdir(path)
    for chr in chr_list:
        if not os.path.exists(os.path.join(path, chr)):
            os.mkdir(os.path.join(path, chr))
    if os.path.exists(os.path.join(path, header)):
        os.remove(os.path.join(path, header))
    for line in vcffile:
        # print(line)
        change_info(line, header, path, chr_list)
    vcffile.close()


if __name__ == '__main__':
    if args.onlychangeInfo:
        # chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
        chr_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
        only_change(args.vcffile, args.onlychangeInfo, chr_list)
    elif args.outputpath:
        sv_class(args.vcffile, path=args.outputpath)
    else:
        sv_class(args.vcffile)

