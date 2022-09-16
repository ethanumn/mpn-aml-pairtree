import argparse
import json
import sys, os
import pandas as pd

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from ssm_columns import *
from mpn_aml_columns import *


def load_ssm(ssm_fn):
    """
    Return an ssm files dataframe (assumes tab delimited)
    """
    return pd.read_csv(ssm_fn, sep="\t")


def load_xls(xls_fn):
    """
    Return an xls file dataframe (assumes only one sheet)
    """
    return pd.read_excel(xls_fn)


def create_subpop_file(ssm_fn, params_fn, xls_fn, xls_out):

    ssm_df = load_ssm(ssm_fn)
    xls_df = load_xls(xls_fn)

    clusters = None

    pop_df_list = []

    with open(params_fn, "r") as params_json:
        params_data = json.load(params_json)
        clusters = params_data["clusters"]

    # iterate through populations in clusters
    for pop in clusters:

        # pull information about clusters from ssm file
        pop_vars = ssm_df.loc[ssm_df[COL_ID].isin(pop)]

        gene_position = pop_vars[COL_NAME].str.split("_", expand=True)

        # obtain gene and position information for each variant in population
        gene_arr, position_arr = gene_position[0], gene_position[1]

        # create new pop df
        pop_df = pd.DataFrame()

        # grab all rows in input xls that are related
        for gene, pos in zip(gene_arr, position_arr):

            pop_df = pop_df.append(xls_df.loc[(xls_df[GENE]==gene) & (xls_df[POSITION]==int(pos))])
            print(pop_df)

        pop_df_list.append(pop_df)

    with pd.ExcelWriter(xls_out) as writer:
        pop_num = 1
        for pop_df in pop_df_list:
            pop_df.to_excel(writer, sheet_name="Pop"+str(pop_num))
            pop_num += 1





def main():

    parser = argparse.ArgumentParser(

        description='Process a custom file type into an .ssm file readable by Pairtree',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-s', '--ssm-fn', help="ssm file to reference")
    parser.add_argument('-p', '--params-fn', help='params files to obtain populations from')
    parser.add_argument('-x', '--xls-fn', help='excel file to pull rows from')
    parser.add_argument('-o', '--out-fn', help='excel file to write out subpopulations to')

    args = parser.parse_args()

    create_subpop_file(args.ssm_fn, args.params_fn, args.xls_fn, args.out_fn)

if __name__ == '__main__':
  main()
