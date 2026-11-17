import pandas as pd
from collections import defaultdict
import sys


def devsion10(x, count=0):
    if x > 100:
        count += 1
        x = x // 10
        return devsion10(x, count)
    else:
        return x , count

def custom_round(number):
    if number < 100000:
        # get the integer and decimal parts of the number
        integer_part , count = devsion10(number)
        decimal_part = integer_part % 10
        # print(integer_part, decimal_part, count)
        integer_part = integer_part // 10
        count += 1
        # round the decimal part to the nearest 10
        if decimal_part < 5:
            return integer_part * 10 ** (count)
        else:
            return (integer_part + 1) * 10** (count)

def sv_info2dict(sv_info_file, min_sv_len=1, max_sv_len=100000):
    # only consider SVs with length > 10bp and < 1Mb
    # sv_sample_dict = defaultdict()
    sv_info_dict = defaultdict(lambda: defaultdict(dict))
    with open(sv_info_file, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('tools'):
                tool_name, sampleID, pos, sv_type, sv_len = line.strip().split('\t')
                if max_sv_len > int(sv_len) > min_sv_len:
                    value = custom_round(float(sv_len))
                    sv_info_dict[sampleID][sv_type][str(value)] = sv_info_dict[sampleID][sv_type].get(str(value), 0) + 1
                    
    return sv_info_dict

def sys_info():
    
    if len(sys.argv)!= 3:
        print('Usage: python sv_len_stats_plot.py sv_info_file output_file')
        sys.exit(1)

if __name__ == '__main__':
    sys_info()
    sv_info_file = sys.argv[1]
    prefix_output_file = sys.argv[2]
    sv_info_dict = sv_info2dict(sv_info_file)
    with open(prefix_output_file+'_count_plot.txt', 'w') as f:
        print('SampleID\tSV_type\tLength\tCount', file=f)
        for sampleID in sv_info_dict:
            for svtype in sv_info_dict[sampleID]:
                for length, count in sv_info_dict[sampleID][svtype].items():
                    print(sampleID, svtype, length, count, sep='\t', file=f)
    
    # df1 = pd.DataFrame(sv_info_dict).reset_index().rename(columns={'index': 'Length'})
    # df1.fillna(0, inplace=True)
    # for col in df1.columns[1:]:
    #     df1[col] = df1[col].apply(lambda x: x / df1[col].sum())

    # df2 =pd.melt(df1, id_vars=['Length'], var_name='sv_type', value_name='Count')
    # df2.to_csv(prefix_output_file+'_density_plot.txt', sep='\t', index=False)