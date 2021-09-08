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


def overwrite_ids(dataframe):
    """
    Overwrite id column after some change has been made which modifies the original order of the dataframe
    """

    # overwrite id column such that matched columns are at the end
    dataframe[COL_ID] = ["s" + str(number) for number in list(range(0, len(dataframe)))]

    return dataframe


def remove_vars_by_vaf(dataframe, op, vaf):
    """
    Breaks apart string representation of variant and total reads,
    and uses these arrays to find which rows (locus with a variant)
    meets the criteria defined by the passed operator and VAF
    """

    rows_to_keep = np.invert([
        any([
        op(int(vr)/int(tr), vaf) for vr, tr in zip(var_reads, total_reads)
        ]) \
        for var_reads, total_reads in zip(dataframe[COL_VAR_READS].apply(lambda row: row.split(",")), dataframe[COL_TOTAL_READS].apply(lambda row: row.split(",")))
    ])



    return overwrite_ids(dataframe[rows_to_keep].copy())


def organize_vars_by_vaf(dataframe, op, vaf, var_read_prob=None):
    """
    Breaks apart string representation of variant and total reads,
    and uses these arrays to find which rows (locus with a variant)
    meets the criteria defined by the passed operator and VAF, then
    uses this to reorganize order of rows by placing those which met the
    criteria at the end. Optional to overwrite var_read_prob for matching group as well.
    """
    rows_at_end = [
        any([
        op(int(vr)/int(tr), vaf) for vr, tr in zip(var_reads, total_reads)
        ]) \
        for var_reads, total_reads in zip(dataframe[COL_VAR_READS].apply(lambda row: row.split(",")), dataframe[COL_TOTAL_READS].apply(lambda row: row.split(",")))
    ]

    # boolean series of non-matched rows
    rows_at_beginning = np.invert(rows_at_end)

    # overwrite var_read_prob for matched rows if it's passed in
    if var_read_prob:
        dataframe.loc[rows_at_end, COL_VAR_READ_PROB] = \
            dataframe.loc[rows_at_end, COL_VAR_READS].apply(lambda row: ", ".join(len(row.split(","))*[var_read_prob]))

    # append the two dataframes together
    dataframe = dataframe[rows_at_beginning].append(dataframe[rows_at_end]).reset_index(drop=True)

    # rewrite ids
    dataframe = overwrite_ids(dataframe)

    return dataframe


def scale_counts(dataframe, *cell_counts):
    """
    Scales var_reads and total_reads by a max_total_reads on a per variant/sample basis,
    such that all((total_reads <= max_total_reads) and (var_reads <= max_total_reads))
    """

    import re, json

    # handle whether or not we've passed in a general cell count estimate or a per sample cell count estimate
    estimated_coverage = None

    if len(cell_counts) == 2:

        if re.search(".*\.params\.json", cell_counts[0]) and re.search(".*\.csv", cell_counts[1]): # if we're reading from a file

            sample_order = json.load(open(cell_counts[0]))["samples"]
            cell_counts_df = load_csv(cell_counts[1])
            estimated_coverage = np.array([cell_counts_df.loc[cell_counts_df["sample"] == sample, "cells"].item() for sample in sample_order]).astype('int64') * 2

        else:
             raise TypeError("Incorrect parameters passed through arguments - expected [.*\.params\.json, .*\.csv]")
    else:
        estimated_coverage = np.array(cell_counts*len(dataframe.iloc[0][COL_TOTAL_READS].split(","))).astype('int64') * 2

    # iterate through all rows and scale var_reads and total_reads using estimated_coverage vector
    for row_idx in range(0, len(dataframe)):

        row = dataframe.iloc[row_idx]

        var_reads_scaled, total_reads_scaled = [], []

        # iterate through all var_read/total_reads per row, scale if necessary, and add values to new vector of scaled values
        for var_reads, total_reads, max_read_count in zip([int(cnt) for cnt in row[COL_VAR_READS].split(",")],
                                                 [int(cnt) for cnt in row[COL_TOTAL_READS].split(",")],
                                                 estimated_coverage
                                                 ):

            factor = 1

            if total_reads > max_read_count:

                factor = total_reads / max_read_count

            var_reads_scaled.append(str(round(var_reads / factor)))
            total_reads_scaled.append(str(round(total_reads / factor)))

        dataframe.loc[row_idx, [COL_VAR_READS, COL_TOTAL_READS]] = [", ".join(var_reads_scaled), ", ".join(total_reads_scaled)]

    return dataframe


def keep_vars_by_name(dataframe, names):
    """
    Returns a dataframe that has only the names given
    """

    return overwrite_ids(dataframe[dataframe[COL_NAME].isin(names)])



def separate_garbage(dataframe, *params_file):
    """
    Returns a dataframe with only the garbage mutations per the given params file
    """
    if re.search(".*\.params\.json", params_file[0]): # if we're reading from a file

        garbage = json.load(open(params_file[0]))["garbage"]

        return dataframe.loc[dataframe[COL_ID] == garbage].copy()

    else:

        return pd.Dataframe()



def pyclone_vi_fmt(dataframe, params):
    """
    Processes ssm dataframe into a tsv that can be used by PyClone-VI
    """

    import json

    samples = json.load(open(params))["samples"]

    pyclone_vi_columns = ["mutation_id", "sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn"]

    pyclone_df = pd.DataFrame(columns=pyclone_vi_columns)

    for row_idx in range(0, len(dataframe)):

        row = dataframe.iloc[row_idx]

        var_reads_scaled, total_reads_scaled = [], []

        # iterate through all var_read, ref_reads, var_read_prob per row
        for var_reads, ref_reads, var_read_prob, sample in zip([int(cnt) for cnt in row[COL_VAR_READS].split(",")],
                                                               [int(cnt) for cnt in row[COL_TOTAL_READS].split(",")],
                                                               [float(prob) for prob in row[COL_VAR_READ_PROB].split(",")],
                                                               samples):
            copy_number = 2 if var_read_prob < 1 else 1
            pyclone_values = [row[COL_NAME], sample, ref_reads, var_reads, copy_number, copy_number, copy_number]
            pyclone_df = pyclone_df.append(dict(zip(pyclone_vi_columns, pyclone_values)), ignore_index=True)

    return pyclone_df
