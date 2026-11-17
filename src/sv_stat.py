# Unix
import pandas as pd
from collections import defaultdict
import os
import argparse
import gzip
import glob
from multiprocessing import Pool

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

def four_sv_stat(filepaht, SVTYPE):
    sv_dict = defaultdict(int)
    f = open_vcf_file(filepaht)
    for line in f:
        if not line.startswith('#'):
            for sv_type in SVTYPE:
                if sv_type in line:
                    sv_dict[sv_type] = sv_dict.get(sv_type, 0) + 1
    for sv_type in SVTYPE:
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
            
        vcf_list = glob.glob(os.path.join(path, "**", sampleID, "**", vcf), recursive=True)
        # print(vcf_list)
        # print(tools_name)
        if len(vcf_list) == 1:
            # print(vcf_list)
            sample_dict = four_sv_stat(vcf_list[0], SVTYPE)
            sample_dict['key1'] = sampleID
            sample_dict['key2'] = tools_name
            tools_list.append(sample_dict)
        elif len(vcf_list) != 1:
            print(f"There are {len(vcf_list)} {tools_name} VCF files for {sampleID}")
    return tools_list

def sample_stat (path, sampleID_list, five_vcf_list, SVTYPE= ["DEL", "INS", "INV", "DUP"]):
    sample_list = []
    # pool = Pool(processes=24)
    for sampleID in sampleID_list:
        # tools_list = pool.apply_async(tools_stat, args=(path, sampleID, five_vcf_list, SVTYPE))
        tools_list = tools_stat(path, sampleID, five_vcf_list, SVTYPE)
        sample_list.extend(tools_list)
    # pool.close()
    # pool.join()
    return sample_list

def list2df(sample_list, prefix, SVTYPE= ["DEL", "INS", "INV", "DUP"]):
    # import openpyxl
    df = pd.DataFrame(sample_list)
    # 将key2作为列名
    df = df.pivot(index='key1', columns='key2', values= SVTYPE)
    # 重置列名
    df.columns = pd.MultiIndex.from_tuples([(j, i) for i, j in df.columns])
    # df = df.swaplevel(axis=1)
    df.sort_index(axis=1, level=0, inplace=True)
    df.dropna(inplace=True, axis=1)
    df.index.name = None
    df.to_excel(prefix + "_sv_stat.xlsx")

#  = ["lumpy/*.genotyped.vcf", "breakdancer/*cfg.SV.output", 'manta/*.manta.sv.vcf','cnvnator/*cnvnator.vcf','metasv/*.SV.vcf.gz']
if __name__ == '__main__':
    args = parser.parse_args()
    path = args.path
    all_sample_list = args.all_sample_list
    # we can add and remove sv type
    SVTYPE = ["DEL", "INS", "INV", "DUP",  "ITX", "CTX"]
    sampleID_list = file2list(all_sample_list)
    
    # statistcs for raw SVs identified by different tools for all samples
    raw_five_vcf_list = ["lumpy/*.genotyped.vcf", "breakdancer/*cfg.SV.output", 'manta/*.manta.sv.vcf','cnvnator/*cnvnator.vcf','metasv/*.SV.vcf.gz']
    raw_sample_sv_list = sample_stat(path, sampleID_list, raw_five_vcf_list, SVTYPE)
    list2df(raw_sample_sv_list, "raw", SVTYPE)
    
    # statistcs for SVs identified by different tools after filtering for all samples
    filter_five_vcf_list = ["metasv/*.lumpy.gt.vcf", "metasv/*manta.gt.vcf", 'metasv/breakdancer.vcf.gz','metasv/cnvnator.vcf.gz','metasv/*.SV.pass.vcf']
    filter_sample_sv_list = sample_stat(path, sampleID_list, filter_five_vcf_list, SVTYPE)
    list2df(filter_sample_sv_list, "filtered", SVTYPE)
