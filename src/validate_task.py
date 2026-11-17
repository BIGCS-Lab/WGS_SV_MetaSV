#coding=utf-8
from __future__ import print_function
import logging
import os
import glob
import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--version', action='version', version='%(prog)s 20240627')
parser.add_argument('--path','-p',type=str,help="the path of the log files",required= True,metavar='')
parser.add_argument('--completed_tasks','-c',type=str,help="the completed tasks file",required= True,metavar='')
parser.add_argument('--err_log','-e',type=str,help="the error log file",required= True,metavar='')
parser.add_argument('--delete_snakemake_dir', action='store_true', default=False, help='remove the .snakemake directory [default: False]')
parser.add_argument('--delete_data', action='store_true', default=False, help='remove the data  or manta directory [default: False]')
args = parser.parse_args()

# logging.basicConfig(level=logging.INFO)
# logging = logging.getlogging()

def validate_done_file(directory_path, file_name):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):
        return True
    else:
        return False

def remaining_task(all_tasks, done_tasks):
    with open(all_tasks, 'r') as f:
        tasks = f.readlines()
    with open(done_tasks, 'r') as f:
        done = f.readlines()
    remaining = [t for t in tasks if t not in done]
    return remaining

def time_difference(time_str1, time_str2):
    from datetime import datetime
    # 定义时间格式
    time_format = "%a %b %d %H:%M:%S CST %Y"
    # 将时间字符串转换为datetime对象
    time1 = datetime.strptime(time_str1, time_format)
    time2 = datetime.strptime(time_str2, time_format)
    # 计算两个时间之间的差异
    time_difference = time2 - time1
    # 将时间差异转换为小时数
    hours_difference = time_difference.total_seconds() // 3600
    return hours_difference

def delete_file(file_list):
    "delete file list"
    for file_path in file_list:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"warning: you have rm the {file_path} file.")
            
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

def rename_file(old_file_path, new_file_path):
    if os.path.exists(new_file_path):
        pass
    else:
        try:
            os.rename(old_file_path, new_file_path)
        except Exception as e:
            pass

def rm_undone_Data_Manta_dir(path, sampleid):
    breakdancer_dir = os.path.join(path, "breakdancer")
    cnvnator_dir = os.path.join(path, "cnvnator")
    lumpy_dir = os.path.join(path, "lumpy")
    manta_dir = os.path.join(path, "manta")

    breakdancer_file = sampleid + ".cfg.SV.output"
    cnvnator_file = sampleid + ".cnvnator.vcf"
    lumpy_file = sampleid + ".genotyped.vcf"
    manta_file = sampleid + ".manta.vcf.gz"
    if (validate_done_file(breakdancer_dir, breakdancer_file) and 
        validate_done_file(cnvnator_dir, cnvnator_file) and 
        validate_done_file(lumpy_dir, lumpy_file) and 
        validate_done_file(manta_dir, manta_file)):
        delete_directory(os.path.join(path, 'data'))
    if validate_done_file(manta_dir, manta_file):
        pass
    else:
        delete_directory(os.path.join(path, 'manta'))
        
    lumpy_splitters_unsortbam = sampleid + ".splitters.unsorted.bam"
    lumpy_splitters_bam = sampleid + ".splitters.bam"
    lumpy_discordants_unsortbam = sampleid + ".discordants.unsorted.bam"
    lumpy_discordants_bam = sampleid + ".discordants.bam"
    lumpy_vcf = sampleid + ".lumpy.vcf"
    if not validate_done_file(lumpy_dir, lumpy_vcf):
        delete_directory(lumpy_dir)
        
    breaddancer_res = sampleid + ".cfg"
    if not validate_done_file(breakdancer_dir, breakdancer_file):
        delete_directory(breakdancer_dir)
        
    if not validate_done_file(cnvnator_dir, cnvnator_file):
        delete_directory(cnvnator_dir)
    
    
    
        
def check_done_tasks(log_file_path, err_log, completed_tasks, delete_snakemake_dir=False, delete_data=False):
    err_log_path = open(err_log, 'a')
    completed_tasks_path = open(completed_tasks, 'a')
    print(file=err_log_path)
    print(f"Checking done tasks in {log_file_path}", file=err_log_path)
    with open(log_file_path, 'r') as f:
        next(f)
        node = next(f).strip()
        task_name = next(f).strip().split(":")[0]
        # done_list = []
        for line in f:
            if line.strip().startswith('Command') and 'completed' in line:
                cmd = line.strip().split(" ")[2]
                path = os.path.dirname(cmd)
                sampleid = path.split("/")[-1]
                file_name = path.split("/")[-1]+"_have_done.txt"
                log_dir = os.path.dirname(log_file_path)
                # check if the "done" file exists in the directory
                if validate_done_file(path, file_name):
                    # logging.info(f"{file_name} exists in {path}, skip.")
                    # done_list.append(path.split("/")[-1])
                    print(f'{log_dir.split("/")[-1]}: {path}', file=completed_tasks_path)
                    old_snakemake_err_txt = os.path.join(path, 'snakemake.err.txt')
                    new_snakemake_err_txt = os.path.join(path, 'snakemake.err.1.txt')
                    if os.path.exists(new_snakemake_err_txt):
                        pass
                    else:
                        rename_file(old_snakemake_err_txt, new_snakemake_err_txt)
                else:
                    print(f"{task_name} has completed in {node}, but {file_name} does not exist in {path}, pls check the following directory.", file=err_log_path)
                    print(f"cd {path}", file=err_log_path)
                    _need_delete_dir = os.path.join(path, '.snakemake')
                    if delete_snakemake_dir:
                        delete_directory(_need_delete_dir)
                    if delete_data:
                        rm_undone_Data_Manta_dir(path, sampleid)
                    # output the errored tasks to the error task file
                    if "30x" in log_file_path:
                        err_task_file = "30x_err_task.sh"
                        with open(err_task_file, 'a') as f:
                            print(f'sh {cmd}', file=f)
                    elif "7x" in log_file_path:
                        err_task_file = "7x_err_task.sh"
                        with open(err_task_file, 'a') as f:
                            print(f'sh {cmd}', file=f)
                    else:
                        err_task_file = "all_err_task.sh"
                        with open(err_task_file, 'a') as f:
                            print(f'sh {cmd}', file=f)

            if line.startswith('end'):
                time1 = line.split('=')[1].strip()
            elif line.strip().startswith('start'):
                time2 = line.split('=')[1].strip()
                
        if "30x" in log_file_path:
            hours_difference = time_difference(time1, time2)
            if hours_difference > 8:
                print(f"The {task_name} in {log_file_path} is running in {node} for more than 7 hours, pls check the following directory.", file=err_log_path)
                print(f"cd {path}", file=err_log_path)
        elif "7x" in log_file_path:
            hours_difference = time_difference(time1, time2)
            if hours_difference > 3:
                print(f"The {task_name} in {log_file_path} is running in {node} for more than 3 hours, pls check the following directory.", file=err_log_path)
                print(f"cd {path}", file=err_log_path)
    err_log_path.close()
    completed_tasks_path.close()
    
if __name__ == '__main__':
    print('Warning: you need to whether to quit the jobs of the remaining tasks in the error log file, because 30x_err_task.sh and 7x_err_task.sh will be rewrited')
    log_list = glob.glob(os.path.join(args.path, '**', 'slurm*.out'), recursive=True)
    delete_file([args.err_log, args.completed_tasks, "30x_err_task.sh", "7x_err_task.sh", "all_err_task.sh"])
    if args.delete_snakemake_dir:
        for log_file in log_list:
            check_done_tasks(log_file, args.err_log, args.completed_tasks, delete_snakemake_dir=True)
    elif args.delete_data:
        for log_file in log_list:
            check_done_tasks(log_file, args.err_log, args.completed_tasks, delete_data=True)
    
    elif args.delete_snakemake_dir and args.delete_data:
        for log_file in log_list:
            check_done_tasks(log_file, args.err_log, args.completed_tasks, delete_snakemake_dir=True, delete_data=True)
    else:
        for log_file in log_list:
            check_done_tasks(log_file, args.err_log, args.completed_tasks)
        
