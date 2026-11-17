#coding=utf-8
from __future__ import print_function
import logging
import os
import glob
import argparse
example_text = '''example:
    python this.py -p path -a all_sample.txt
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--version', action='version', version='%(prog)s the version is 20240704')
parser.add_argument('--path','-p',type=str,help="the path of the output files",required= True,metavar='')
parser.add_argument('--all_sample', '-a', type=str,help="the all sample id",required= True,metavar='')
parser.add_argument('--delete', action='store_true', default=False, help='remove the data  or manta directory [default: False]')
args = parser.parse_args()


import shutil

def delete_directory(directory_path):
    """
    delete directory recursively
    """
    try:
        shutil.rmtree(directory_path)
        print(f'warning: you have rm the {directory_path} directory.')
    except Exception as e:
        pass

def delete_file(file_list):
    "delete file list"
    for file_path in file_list:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f'warning: you have rm the {file_path} file.')


if __name__ == '__main__':
    path = args.path
    all_sample = args.all_sample
    sample_list = []
    with open(all_sample, 'r') as f:
        for line in f:
            sample_list.append(line.strip().split('\t')[0])
    
    done_task = open('done_task.txt', 'w')
    undone_task = open('undone_task.txt', 'w')

    for sampleid in sample_list:
        done_file = os.path.join(path, sampleid, sampleid+"_have_done.txt")
        breakdancer_res = os.path.join(path, sampleid, "breakdancer",  sampleid+ ".cfg.SV.output")
        lumpy_res = os.path.join(path,sampleid,  "lumpy",  sampleid+ ".genotyped.vcf")
        manta_res = os.path.join(path,sampleid,  "manta",  sampleid+ ".manta.sv.vcf")
        manta_gt_vcf = os.path.join(path, sampleid, "metasv",  sampleid+ ".manta.gt.vcf")
        lumpy_gt_vcf = os.path.join(path,sampleid,  "metasv",  sampleid+ ".lumpy.gt.vcf")
        metasv_res = os.path.join(path, sampleid, "metasv",  sampleid+ ".SV.pass.vcf")
        bam_index = os.path.join(path, sampleid, "data",  sampleid+ ".bam.bai")
        bam = os.path.join(path, sampleid, "data",  sampleid+ ".bam")
        if (os.path.exists(breakdancer_res) and os.path.exists(lumpy_res) 
            and os.path.exists(manta_res) and os.path.exists(manta_gt_vcf) 
            and os.path.exists(lumpy_gt_vcf) and os.path.exists(metasv_res) 
            and os.path.exists(done_file)):
            done_task.write(sampleid+'\t'+os.path.join(path, sampleid)+'\n')
        else:
            undone_task.write(sampleid+'\t'+os.path.join(path, sampleid)+'\n')
            if args.delete:
                if os.path.exists(bam_index) and os.path.exists(bam):
                    pass
                else:
                    delete_file([bam_index, bam])
                    
                if os.path.exists(manta_res):
                    delete_directory(os.path.join(path, sampleid, "manta" ))

    done_task.close()
    undone_task.close()

