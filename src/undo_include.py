#coding=utf-8
from __future__ import print_function 
from __future__ import division
import argparse
import os
import glob
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--all_tasks','-a',type=str,help="",required= True,metavar='')
# parser.add_argument('--done_tasks', '-d',type= str,help="",required= True,metavar='')
parser.add_argument('--path', '-p',type= str,help="",required= True,metavar='')
# parser.add_argument('--prefix', '-p',type= str,help="",required= True,metavar='')
args = parser.parse_args()

# def done_tasks(done_file):
#     with open(done_file, 'r') as done_file:
#         done_list = [line.strip() for line in done_file.readlines()]
#     return done_list

# def get_done_tasks(path, suffix="*have_done.txt"):
#     done_list = []
#     done_files = glob.glob(os.path.join(path, '**', suffix), recursive=True)
#     count = len(done_files)
#     print("Find {} done files".format(count))
#     if count != 0:
#         for done_file in done_files:
#             with open(done_file, 'r') as done_file:
#                 done_list.extend([line.strip().split(' ')[0] for line in done_file.readlines()])
#     return done_list

def get_done_tasks(path, all_info, suffix="*have_done.txt"):
    done_list = []

def delete_file(file_list):
    "delete file list"
    for file_path in file_list:
        if os.path.exists(file_path):
            os.remove(file_path)
            
if __name__ == '__main__':
    done_list = get_done_tasks(args.path)
    delete_file(['7X_todo.sh', '30X_todo.sh'])
    # with open(args.prefix + '_todo.sh', 'w') as todo_file:
    with open(args.all_tasks, 'r') as all_file:
        for line in all_file:
            lines = line.strip().split('\t')
            sampleid = lines[0].strip()
            depth = lines[1].strip()
            if sampleid not in done_list:
                with open(depth + '_todo.sh', 'a') as todo_file:
                    todo_file.write(line)
