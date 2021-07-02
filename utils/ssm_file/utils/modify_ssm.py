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


def scale_counts(dataframe, max_total_reads):
    """
    Scales var_reads and total_reads by a max_total_reads on a per variant/sample basis,
    such that all((total_reads <= max_total_reads) and (var_reads <= max_total_reads))
    """

    for row_idx in range(0, len(dataframe)):

        row = dataframe.iloc[row_idx]

        var_reads_scaled, total_reads_scaled = [], []

        # iterate through all var_read/total_reads per row, scale if necessary, and add values to new vector of scaled values
        for var_reads, total_reads in zip([int(cnt) for cnt in row[COL_VAR_READS].split(",")], [int(cnt) for cnt in row[COL_TOTAL_READS].split(",")]):

            factor = 1

            if total_reads > max_total_reads:

                factor = total_reads / max_total_reads

            var_reads_scaled.append(str(round(var_reads / factor)))
            total_reads_scaled.append(str(round(total_reads / factor)))

        dataframe.loc[row_idx, [COL_VAR_READS, COL_TOTAL_READS]] = [", ".join(var_reads_scaled), ", ".join(total_reads_scaled)]

    return dataframe
