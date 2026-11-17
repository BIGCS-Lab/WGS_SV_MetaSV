from collections import defaultdict
import os
import re
from multiprocessing import Pool
import argparse
import pandas as pd
import numpy as np

def parse_line(line):
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE = re.compile(r'\bSVTYPE=(-?\w+)\b')
    chrom, pos, id, ref, alt, qual, filter, info, *gtlist = line.strip().split('\t')
    # if chrom in chr_list:
    try:
        end = ENDpos.findall(info)[0]
    except:
        pass
    try:
        svlen = SVLEN.findall(info)[0]
        end = ENDpos.findall(info)[0]
    except:
        svlen = int(end) - int(pos)
    svtype = SVTYPE.findall(info)[0]
    loci = str(chrom) + ':' + str(pos) + '-' + str(end)
    return svtype, loci, abs(int(svlen))
            # print(loci, svtype, svlen, file=headerf)

def vcff2dict(vcffile):
    svdict = defaultdict(lambda: defaultdict(int))
    tmp_list= []
    with open(vcffile, 'r') as f:
        pool = Pool(processes=24)
        for line in f:
            if line.startswith('#'):
                pass
            else:
                _list = pool.apply_async(parse_line, args=(line,))
                tmp_list.append(_list)
        for _list in tmp_list:
            svtype, loci, svlen = _list.get()
            svdict[svtype][loci]=svlen
        pool.close()
        pool.join()
    return svdict

def main():
    parser = argparse.ArgumentParser(description='sv_stats_individual.py')
    parser.add_argument('-i', '--input', required=True, help='input vcf file')
    parser.add_argument('-p', '--prefix', required=True, help='prefix of output files')
    args = parser.parse_args()
    inputf = args.input
    prefix = args.prefix
    svdict = vcff2dict(inputf)
    with open(prefix+'_sv_stats.txt', 'w') as f:
        print('#svtype\tloci\tsvlen', file=f)
        for svtype in svdict:
            for loci in svdict[svtype]:
                svlen = svdict[svtype][loci]
                f.write(svtype + '\t' + loci + '\t' + str(svlen) + '\n')
    df_sv = pd.read_csv(prefix+'_sv_stats.txt', sep='\t', header=0)
    df_sv.groupby(['#svtype']).size().reset_index(name='count').to_csv(prefix+'_sv_type_count.txt', sep='\t', index=False)
    df_len_stats = pd.pivot_table(df_sv, index=['#svtype'], values=['svlen'], aggfunc=[np.sum, np.max, np.min, np.mean, np.median])
    df_len_stats.to_csv(prefix+'_sv_len_stats.txt', sep='\t', index=True)

if __name__ == '__main__':
    main()