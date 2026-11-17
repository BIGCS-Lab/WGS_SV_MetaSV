import sys
from pathlib import Path

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
                f.write(f'{task}\n')
    with open(f'{dirname}/work.{num_tasks+1}.sh', 'w') as f:
        f.write(f'#!/bin/bash\n')
        f.write(f'# yhbatch -N 1 -n 24 -p rhenv\n')
        for task in tasks[num_tasks*chunk_size:]:
            f.write(f'{task}\n')

if __name__ == '__main__':
    task_file = sys.argv[1]
    num_tasks = int(sys.argv[2])
    split_task(task_file, num_tasks-1)