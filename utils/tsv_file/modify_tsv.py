import pandas as pd
import numpy as np
import sys, os, math

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from ssm_columns import *


def load_ssm(in_file):
    """
    Return an ssm files dataframe (assumes tab delimited)
    """
    return pd.read_csv(in_file, sep="\t")


def load_csv(csv_file):
    """
    Return a dataframe created from a csv file
    """
    return pd.read_csv(csv_file, sep=",")


def save_ssm(dataframe, out_file):
    """
    Given a dataframe and out_file, saves a tab delimited dataframe
    """
    dataframe.to_csv(out_file, sep="\t", index=False)


def match_tsv_to_ssm(tsv_fn, ssm_fn, params_fn):
    """
    Matches tsv data to ssm data to update clusters in params file
    """
    import json, re

    tsv_df = load_ssm(tsv_fn)
    ssm_df = load_ssm(ssm_fn)
    params_data = None

    clusters = []

    # read through cluster ids in tsv
    for cluster_id in tsv_df["cluster_id"].unique():

        cl = []

        # find all mutations in cluster
        mut_in_cluster = tsv_df.loc[tsv_df["cluster_id"] == cluster_id, "mutation_id"].unique()

        # get the id from the ssm file for each mutation in the cluster
        for mut_id in mut_in_cluster:

            cl.append(ssm_df.loc[ssm_df[COL_NAME] == mut_id, COL_ID].mode()[0])

        clusters.append(cl)


    with open(params_fn, "r") as params_json:
        params_data = json.load(params_json)

    if params_data:
        params_data["clusters"] = clusters

    with open(params_fn, "w") as params_json:
        json_out = json.dumps(params_data)
        params_json.write(json_out)
