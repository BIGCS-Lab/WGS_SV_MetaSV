from __future__ import print_function
import json
import sys
import os
import pandas as pd
import glob

def args_parser():
    if len(sys.argv)!= 3:
        print("Usage: python sum_json.py <json_file_path> <output_file_path>")
        sys.exit(1)

def sum_json(json_file_path, output_file_path):
    json_files = glob.glob(json_file_path + "/*/summary.json")
    for i, json_file in enumerate(json_files):
        if i < 1:
            with open(json_file, "r") as f:
                name = json_file.split("/")[-2]
                data = json.load(f)
                data.popitem()
                df1 =pd.DataFrame.from_dict(data,orient='index')
                df1.rename(columns={0:name},inplace=True)
                # print(df1)
        else:
            with open(json_file, "r") as f:
                data = json.load(f)
                data.popitem()
                name = json_file.split("/")[-2]
                df2 =pd.DataFrame.from_dict(data,orient='index')
                df2.rename(columns={0:name},inplace=True)
                # print(df2)
            df1 = pd.merge(df1, df2, left_index=True, right_index=True)
    df1.to_csv(output_file_path)
            
if __name__ == "__main__":
    args_parser()
    json_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    sum_json(json_file_path, output_file_path)
