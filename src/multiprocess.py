#coding=utf-8
from __future__ import print_function
import argparse
import os
import logging
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--all_cmd','-i',type=str,help="the file containing all the commands to be executed",required= True,metavar='')
parser.add_argument('--Number_of_cores', '-c',type= str,help="how many cores to use at the same time",required= True,metavar='')
parser.add_argument('--Number_of_cmd', '-n',type= str,help="how many commands to execute at the same time",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="the prefix of the cmd file",required= True,metavar='')

args = parser.parse_args()

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
        logger.info(f"File '{file_path}' has been deleted and will be recreated next time.")
    else:
        logger.info(f"File '{file_path}' does not exist and will be recreated.")
            
def split_list(lst, n):
    """
    split the lst into N parts, each part has the same length as the original lst
    """
    k, m = divmod(len(lst), n)
    return [lst[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]

def split_cmd(cmd_file,num_of_cores, num_of_cmd, prefix="work"):
    with open(cmd_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        # split the cmd lines into N parts
        part_list = split_list(lines, num_of_cores)
        for i, cmd_list in enumerate(part_list):
            # print(cmd_list, i)
            delete_file(f"{prefix}.{str(i+1)}.sh")
            with open(f"{prefix}.{str(i+1)}.sh", "w") as f:
                
                print(f"#!/bin/bash\n# yhbatch -N 1 -n 24 -p rhenv", file=f)
                print(f'echo "start = $(date)"', file=f)
                print(f'echo "$(hostname)"', file=f)
                num = 1
                for index in range(0, len(cmd_list), num_of_cmd):
                    # get 5 elements at a time
                    current_chunk = cmd_list[index:index+num_of_cmd]
                    print(f"echo '{prefix}.{str(i+1)}.sh: Processing chunk {num} of {len(range(0, len(cmd_list), num_of_cmd))} ...'\n", file=f)
                    print(f'echo "start = $(date)"', file=f)
                    print(f'echo "$(hostname)"', file=f)
                    # Remind which cmd is about to run
                    for cmd in current_chunk:
                        print(f"echo 'Processing {cmd} is going and pls wait...'", file=f)
                        
                    for cmd in current_chunk:
                        print(f"{cmd} &", file=f)
                    print("wait", file=f)
                    print("sleep 1", file=f)
                    
                    # Remind which cmd was completed
                    for cmd in current_chunk:
                        print(f"echo 'Command {cmd} completed and pls check the output file'", file=f)
                    print(f'echo "end = $(date)"', file=f)
                    print(f"echo 'Finished processing chunk {num} of {len(range(0, len(cmd_list), num_of_cmd))}'\n", file=f)
                    print(f'echo ', file=f)
                    print(f'echo "----------------------------------------------------------------------------------------"', file=f)
                    num += 1
                print(f'echo "end = $(date)"', file=f)

if __name__ == '__main__':
    split_cmd(args.all_cmd,int(args.Number_of_cores),int(args.Number_of_cmd), args.prefix)