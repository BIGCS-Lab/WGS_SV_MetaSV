import pandas as pd
from collections import defaultdict
import os
import argparse
import gzip
import glob
from multiprocessing import Pool
import re

example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--path','-p',type=str,help="the path of the file list",required= True,metavar='')
parser.add_argument('--all_sample_list', '-i',type= str,help="the file list of all samples",required= True,metavar='')

def file2list(file_path):
    sampleID_list = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines = line.strip().split('\t')
                sampleID_list.append(lines[0])
    return sampleID_list

def open_vcf_file(filename):
    
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

# def four_sv_stat(filepaht, SVTYPE):
#     sv_dict = defaultdict(int)
#     f = open_vcf_file(filepaht)
#     for line in f:
#         if not line.startswith('#'):
#             for sv_type in SVTYPE:
#                 if sv_type in line:
#                     sv_dict[sv_type] = sv_dict.get(sv_type, 0) + 1
#     for sv_type in SVTYPE:
#         if sv_type not in sv_dict:
#             sv_dict[sv_type] = None
#     f.close()
#     return sv_dict

def four_sv_stat(filepath, tools_name, SVTYPEs):
    sv_dict = defaultdict(int)
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE = re.compile(r'\bSVTYPE=(-?\w+)\b')
    f = open_vcf_file(filepath)
    if tools_name == "breakdancer" and filepath.endswith(".output"):
        for line in f:
            if not line.startswith('#'):
                for sv_type in SVTYPEs:
                    if sv_type in line.strip().split('\t')[6]:
                        sv_dict[sv_type] = sv_dict.get(sv_type, 0) + 1
    elif tools_name in ["lumpy", "breakdancer", "manta", "cnvnator", "metasv"]:
        for line in f:
            if not line.startswith('#'):
                for sv_type in SVTYPEs:
                    chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                    if sv_type in SVTYPE.findall(info)[0]:
                        sv_dict[sv_type] = sv_dict.get(sv_type, 0) + 1
    for sv_type in SVTYPEs:
        if sv_type not in sv_dict:
            sv_dict[sv_type] = None
    f.close()
    return sv_dict

def tools_stat(path, sampleID, five_vcf_list, SVTYPE, tools_names= ["lumpy", "breakdancer", "manta", "cnvnator", "metasv"]):
    tools_list = []
    for vcf in five_vcf_list:
        # to get the tool name
        # tools_name = os.path.dirname(vcf)
        file_name = os.path.basename(vcf)
        dir_name = os.path.dirname(vcf)
        for tool in tools_names:
            if tool in file_name:
                tools_name = tool
                break
            elif tool in dir_name:
                tools_name = tool
                break
            
        vcf_list = glob.glob(os.path.join(path, sampleID, vcf), recursive=True)
        # print(vcf_list)
        # print(tools_name)
        if len(vcf_list) == 1:
            # print(vcf_list)
            sample_dict = four_sv_stat(vcf_list[0],tools_name, SVTYPE)
            sample_dict['key1'] = sampleID
            sample_dict['key2'] = tools_name
            tools_list.append(sample_dict)
        elif len(vcf_list) != 1:
            print(f"There are {len(vcf_list)} {tools_name} VCF files for {sampleID}")
    return tools_list

def sample_stat (path, sampleID_list, five_vcf_list, SVTYPE= ["DEL", "INS", "INV", "DUP"]):
    sample_list = []
    tmp = []
    pool = Pool(processes=24)
    # n = 1
    total_num = len(sampleID_list)
    for sampleID in sampleID_list:
        # print(f"processing {sampleID} {n}/{total_num}")
        tools_list = pool.apply_async(tools_stat, args=(path, sampleID, five_vcf_list, SVTYPE))
        # sample_list.extend(tools_list.get())
        tmp.append(tools_list)
        # n += 1
    for id in tmp:
        sample_list.extend(id.get())    
    pool.close()
    pool.join()
    tmp = []
    return sample_list

def list2df(sample_list, prefix, SVTYPE= ["DEL", "INS", "INV", "DUP"]):
    # import openpyxl
    df = pd.DataFrame(sample_list)
    long_df2 = pd.melt(df, id_vars=['key1', 'key2'], value_vars=SVTYPE, var_name='SVTYPE', value_name='Count')
    long_df2.dropna(inplace=True)
    # long_df2.to_csv(prefix + "_sv_stat_long.csv", index=False)
    long_df2.to_csv(prefix + "_sv_stat_long.txt", sep='\t', index=False)
    # 将key2作为列名
    df = df.pivot(index='key1', columns='key2', values= SVTYPE)
    # 重置列名
    df.columns = pd.MultiIndex.from_tuples([(j, i) for i, j in df.columns])
    # df = df.swaplevel(axis=1)
    df.sort_index(axis=1, level=0, inplace=True)
    df = df.dropna(axis=1, how='all')
    # df.fillna(0, inplace=True)
    df.index.name = None
    # df.to_csv(prefix + "_sv_stat.csv")
    df.to_csv(prefix + "_sv_stat.txt", sep='\t')

#  = ["lumpy/*.genotyped.vcf", "breakdancer/*cfg.SV.output", 'manta/*.manta.sv.vcf','cnvnator/*cnvnator.vcf','metasv/*.SV.vcf.gz']
if __name__ == '__main__':
    args = parser.parse_args()
    path = args.path
    all_sample_list = args.all_sample_list
    # we can add and remove sv type
    SVTYPE = ["DEL", "INS", "INV", "DUP",  "ITX", "CTX"]
    sampleID_list = file2list(all_sample_list)
    
    # statistcs for raw SVs identified by different tools for all samples
    raw_five_vcf_list = ["lumpy/*.lumpy.vcf", "breakdancer/*cfg.SV.output", 'manta/*.manta.sv.vcf','cnvnator/*cnvnator.vcf','metasv/*.SV.vcf.gz']
    raw_sample_sv_list = sample_stat(path, sampleID_list, raw_five_vcf_list, SVTYPE)
    list2df(raw_sample_sv_list, "raw", SVTYPE)
    print("raw sv statistics done")
    
    # statistcs for SVs identified by different tools after filtering for all samples
    filter_five_vcf_list = ["metasv/lumpy.vcf.gz", "metasv/*manta.gt.vcf", 'metasv/breakdancer.vcf.gz','metasv/cnvnator.vcf.gz','metasv/*.SV.pass.vcf']
    filter_sample_sv_list = sample_stat(path, sampleID_list, filter_five_vcf_list, SVTYPE)
    list2df(filter_sample_sv_list, "filtered", SVTYPE)
    print("filtered sv statistics done")
