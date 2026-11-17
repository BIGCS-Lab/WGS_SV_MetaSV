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
parser.add_argument('--chose', '-c',type= str,help=
                    " [ filtered or raw ] that you need to chose, which means the type of SVs: raw or filtered",
                    required= True,metavar='')

def open_vcf_file(filename):
    '''
    open multiple format of vcf file
    '''
    import gzip
    _vcf_f = None
    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "rt")
    else:
        _vcf_f = open(filename, "r")
    return _vcf_f

def sv_len_stat(file, tools ,sampleID, sv_type = ["DEL", "INS", "INV", "DUP",  "ITX", "CTX"]):
    alist = []
    openf = open_vcf_file(file)
    ENDpos =re.compile(r'\bEND=(-?\d+)\b')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\b')
    SVTYPE = re.compile(r'\bSVTYPE=(-?\w+)\b')
    if tools == "breakdancer" and file.endswith(".output"):
        for line in openf:
            if not line.startswith('#'):
                for sv_type_i in sv_type:
                    if sv_type_i in line.strip().split('\t')[6]:
                        sv_len = abs(int(line.split('\t')[7]))
                        start_pos = str(line.split('\t')[0]) + ":" + str(line.split('\t')[1])
                        # print(tools, sampleID, start_pos, sv_type_i, sv_len, sep='\t', file=output_file)
                        alist.append([tools, sampleID, start_pos, sv_type_i, sv_len])
                        
    elif tools in ["lumpy", "breakdancer", "manta", "cnvnator", "metasv"]:
        for line in openf:
            if not line.startswith('#'):
                for sv_type_i in sv_type:
                    chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
                    if sv_type_i in SVTYPE.findall(info)[0]:
                        try:
                            sv_len = abs(int(SVLEN.findall(info)[0]))
                        except:
                            end = ENDpos.findall(info)[0]
                            sv_len = abs(int(pos) - int(end))
                        if int(sv_len) == 0:
                            end = ENDpos.findall(info)[0]
                            sv_len = abs(int(pos) - int(end))
                        start_pos = chrom + ":" + pos
                        # print(tools, sampleID, start_pos, sv_type_i, sv_len, sep='\t', file=output_file)
                        alist.append([tools, sampleID, start_pos, sv_type_i, sv_len])
    else:
        print("Error: tools not supported")
        # return None
    openf.close()
    return alist

def mul_process(path, five_vcf_list, tools_names, sampleID, sv_type):
    alist = []
    for vcf in five_vcf_list:
        file_name = os.path.basename(vcf)
        dir_name = os.path.dirname(vcf)
        for tool in tools_names:
            if tool in file_name:
                tools_name = tool
                break
            elif tool in dir_name:
                tools_name = tool
        # vcf_file = glob.glob(os.path.join(path, sampleID, tools_name, vcf), recursive=True)
        vcf_file = glob.glob(os.path.join(path, sampleID, vcf), recursive=True)
        if len(vcf_file) == 1:
            tmp_list = sv_len_stat(vcf_file[0], tools_name ,sampleID, sv_type)
            alist.extend(tmp_list)
            
        elif len(vcf_file) != 1:
            print(f"There are {len(vcf_file)} {tools_name} VCF files for {sampleID}")
    return alist

def all2stat(path, five_vcf_list, sampleID_list,  
             tools_names= ["lumpy", "breakdancer", "manta", "cnvnator", "metasv"], 
             sv_type = ["DEL", "INS", "INV", "DUP",  "ITX", "CTX"]):
    # with open(os.path.join(path, outfilename), "w") as output_file:
    all_list = []
    tmp_list = []
    # n = 0
    pool = Pool(processes=24)
    for sampleID in sampleID_list:
        # n += 1
        total_num = len(sampleID_list)
        # print(f"processing {sampleID} {n}/{total_num}")
        _alist = pool.apply_async(mul_process, args=(path, five_vcf_list, tools_names, sampleID, sv_type))
        tmp_list.append(_alist)
    for _list in tmp_list:
        # all_list.extend(_list.get())
        for _item in _list.get():
            print(*_item, sep="\t")
    pool.close()
    pool.join()
    tmp_list = []
    # print(all_list)
    return all_list
            

if __name__ == '__main__':
    args = parser.parse_args()
    path = args.path
    all_sample_list = args.all_sample_list
    # filtered_or_raw = args.chose
    with open(all_sample_list, "r") as f:
        sampleID_list = [line.strip().split()[0] for line in f if not line.startswith("#")]
    print('tools', 'sampleID', 'start_pos', 'sv_type', 'sv_len', sep="\t")
    five_vcf_list = args.chose.split(",")
    all2stat(path, five_vcf_list, sampleID_list)
    
    # if filtered_or_raw == "raw":
        
    #     # statistcs for raw SVs identified by different tools for all samples
    #     raw_five_vcf_list = ["lumpy/*.genotyped.vcf", "breakdancer/*cfg.SV.output", 'manta/*.manta.sv.vcf','cnvnator/*cnvnator.vcf','metasv/*.SV.vcf.gz']
    #     # outfilename = "raw_sv_len_stat.csv"
    #     outfilename = "raw_sv_len_stat.txt"
    #     raw_list = all2stat(path, raw_five_vcf_list, sampleID_list, outfilename)
    #     # with open(outfilename, "w") as output_file:
    #     # # print('tools', 'sampleID', 'start_pos', 'sv_type', 'sv_len', sep=",", file=output_file)
    #     #     print('tools', 'sampleID', 'start_pos', 'sv_type', 'sv_len', sep="\t", file=output_file)  
    #     #     for _list in raw_list:
    #     #         print(*_list, sep="\t", file=output_file)

    # elif filtered_or_raw == 'filtered':
    #     # # statistcs for SVs identified by different tools after filtering for all samples
    #     filter_five_vcf_list = ["metasv/*.lumpy.gt.vcf", "metasv/*manta.gt.vcf", 'metasv/breakdancer.vcf.gz','metasv/cnvnator.vcf.gz','metasv/*.SV.pass.vcf']
    #     outfilename = "filter_sv_len_stat.txt"
    #     filter_list = all2stat(path, filter_five_vcf_list, sampleID_list, outfilename)
    #     # # with open(outfilename, "w") as output_file:
    #     # # # print('tools', 'sampleID', 'start_pos', 'sv_type', 'sv_len', sep=",", file=output_file)
    #     # #     print('tools', 'sampleID', 'start_pos', 'sv_type', 'sv_len', sep="\t", file=output_file)  
    #     # #     for _list in filter_list:
    #     # #         print(*_list, sep="\t", file=output_file)
    