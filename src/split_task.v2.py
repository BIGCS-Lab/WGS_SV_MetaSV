from pathlib import Path
import argparse
example_text = '''example: python * '''
parser = argparse.ArgumentParser(description="The script is .",formatter_class= argparse.RawTextHelpFormatter, usage = '%(prog)s [-h]', epilog=example_text)

parser.add_argument('--task_file','-f',type=str,help="", required= True,metavar='')
parser.add_argument('--num_tasks', '-n',type= int,help="the template of config.yaml",required= True,metavar='')

def split_task(task_file, num_tasks):
    tasks = []
    with open(task_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            tasks.append(line.strip())
    num_tasks = min(num_tasks, len(tasks))
    chunk_size = len(tasks) // num_tasks
    remainder = len(tasks) % num_tasks
    chunks = [tasks[i:i+chunk_size] for i in range(0, len(tasks), chunk_size)]
    dirname = task_file.strip('.sh')
    Path(dirname).mkdir(parents=True, exist_ok=True)
    for i, chunk in enumerate(chunks):
        with open(f'{dirname}/work.{i+1}.sh', 'w') as f:
            f.write(f'#!/bin/bash\n')
            f.write(f'# yhbatch -N 1 -n 24 -p rhenv\n')
            for task in chunk:
                f.write(f'{task} &\n')
            f.write(f'wait\n')
    with open(f'{dirname}/work.{num_tasks+1}.sh', 'w') as f:
        f.write(f'#!/bin/bash\n')
        f.write(f'# yhbatch -N 1 -n 24 -p rhenv\n')
        for task in tasks[num_tasks*chunk_size:]:
            f.write(f'{task} &\n')
        f.write(f'wait\n')

if __name__ == '__main__':
    args = parser.parse_args()
    task_file = args.task_file
    num_tasks = args.num_tasks
    split_task(task_file, num_tasks-1)
