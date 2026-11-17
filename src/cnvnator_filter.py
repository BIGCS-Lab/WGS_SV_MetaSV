import pandas as pd
import numpy as np
import sys

def cnvnator_filter(cnvnator, output_file):
    cnvnator = pd.read_csv(cnvnator, sep='\t', header=None)
    cnvnator.columns = ['svtype', 'region',  'length', 'depth', 'eval_1', 'eval_2','eval_3','eval_5', 'q0']
    #cnvnator_filter = cnvnator[(cnvnator['eval_1']< 0.05) & (cnvnator['eval_2']< 0.05) & (cnvnator['q0']< 0.5) & (cnvnator['q0'] != -1)] 
    cnvnator_filter = cnvnator[(cnvnator['eval_1']< 0.05) & (cnvnator['q0']< 0.5) & (cnvnator['q0'] != -1)] 
    cnvnator_filter.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == '__main__':
    cnvnator_file = sys.argv[1]
    output_file = sys.argv[2]
    cnvnator_filter(cnvnator_file, output_file)
