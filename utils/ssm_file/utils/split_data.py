import pandas as pd
import sys, os
import argparse
import json

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from ssm_columns import *
from mpn_aml_columns_txt import *


def read_params(params_fn):
    """Reads in clusters/garabge from .params.json file"""
    clusters = []
    garbage = []
    
    with open(params_fn, 'r') as params:
        
        data = json.load(params)
        clusters = data[CLUSTERS]
        garbage = data[GARBAGE]

    return clusters, garbage
    

def read_fn(fn):
    """
    Read original data file 
    """
    
    df = pd.DataFrame()

    if fn:
        
        file_ext = fn.split(".")[-1]
            
        if file_ext == "txt" or file_ext == "ssm":
            df = pd.read_csv(fn, sep="\t", header=0)

    return df
    

def main():
    parser = argparse.ArgumentParser(

        description='Break apart original data based on .ssm and .params.json.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-s', '--ssm-files', nargs='+', help='List of ssm files', required=True)
    parser.add_argument('-p', '--param-files', nargs='+', help='List of params.json files', required=True)
    parser.add_argument('-f', '--data-files', nargs='+', help='List of data files', required=True)
    parser.add_argument('-d', '--directories', nargs='+', help='Directory to read/write files from.')


    args = parser.parse_args()
    assert (len(args.ssm_files) == len(args.param_files)) and (len(args.ssm_files) == len(args.data_files)), "Number of ssm, params, and data files must match"
    assert len(args.ssm_files) == len(args.directories) or len(args.directories) == 1, "Number of directories must match number of other input files, or be a single directory"
           
           
    for i, (ssm_fn, params_fn, data_fn) in enumerate(zip(args.ssm_files, args.param_files, args.data_files)):
        
        if len(args.directories) == 1:
            ssm_df = read_fn(args.directories[0] + ssm_fn)
            clusters, garbage = read_params(args.directories[0] + params_fn)
            data_df = read_fn(args.directories[0] + data_fn)
            
        else:
            ssm_df = read_fn(args.directories[i] + ssm_fn)
            clusters, garbage = read_params(args.directories[i] + params_fn)
            data_df = read_fn(args.directories[i] + data_fn)
            
            
        dir = data_fn.split(".")[0]
        os.mkdir(dir)

        # extract clusters data
        for i, c in enumerate(clusters):
                    
            df = pd.DataFrame()

            for id in c:
                chr, pos = ssm_df.loc[ssm_df[COL_ID] == id][COL_NAME].values[0].split("_")
                # print(data_df.loc[(data_df[CHR] == chr) & (data_df[START] == pos)])
                df = df.append(data_df.loc[(data_df[CHR] == chr) & (data_df[START] == int(pos))])
            
            df.to_csv(dir + "/" + "cluster%d.txt" % (i+1), sep="\t", index=False)
                
        
        # extract garbage data
        df = pd.DataFrame()

        for id in garbage:
            chr, pos = ssm_df.loc[ssm_df[COL_ID] == id][COL_NAME].values[0].split("_")
            # print(data_df.loc[(data_df[CHR] == chr) & (data_df[START] == pos)])
            df = df.append(data_df.loc[(data_df[CHR] == chr) & (data_df[START] == int(pos))])
        
        df.to_csv(dir + "/" + "garbage.txt", sep="\t", index=False)   





if __name__ == '__main__':
  main()