import unittest
import os, sys
import argparse

import numpy as np
import pandas as pd


sys.path.append(os.environ["UTILS_DIR"] + "/common")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'xls_aggregators'))

from mpn_aml_columns import *
from mpn_aml_aggregator import MPN_AML_Aggregator


class MPN_AML_Processor_Tests(unittest.TestCase):
    """
    Test cases for MPN_AML_Aggregator for aggregating xls spread sheets for mpn_aml project.
    use 'python3 test_mpn_aml_aggregator.py' to run the test suite
    """
    def setUp(self):

        # set up using example in-file
        self.primary_xls = [os.environ["DATA_DIR"] + "/example/" + "example.primary.xlsx", "Sheet1"]
        self.calls_xls = [os.environ["DATA_DIR"] + "/example/" + "example.calls.xlsx", "Sheet1"]
        self.populations_xls = [os.environ["DATA_DIR"] + "/example/" + "example.populations.xlsx", "Sheet1", None]

        # set up for mpn-aml example
        # self.primary_xls = [os.environ["DATA_DIR"] + "/mpn_aml_test/" + "mats08_primary.xlsx", "SNVs"]
        # self.calls_xls = [os.environ["DATA_DIR"] + "/mpn_aml_test/" + "mats08_calls.xlsx", "Sheet1"]
        # self.populations_xls = [os.environ["DATA_DIR"] + "/mpn_aml_test/" + "mats08_primary.xlsx", "Populations", None]

        self.test_aggregator = MPN_AML_Aggregator(primary_xls = self.primary_xls,
                                                calls_xls = self.calls_xls,
                                                populations_xls = self.populations_xls,
                                                aggregated_xls = [],
                                                write_xls_file = False)

        self.on_list = [CHR_POS, CHR, POSITION, SAMPLE_NAMES, ALT_DEPTH, REF_DEPTH]
        self.on_list_no_alt = [CHR_POS, CHR, POSITION, SAMPLE_NAMES, REF_DEPTH]

        self.primary_df = self.test_aggregator.primary_df
        self.calls_df = self.test_aggregator.calls_df
        self.aggregated_df = self.test_aggregator.aggregated_df

        # inner join between primary_df and aggregated_df
        self.shared_rows_aggregated_primary = pd.merge(self.primary_df, self.aggregated_df, how="inner", on=self.on_list, suffixes=("_x", ""))


        # inner join between aggregated_df and calls_df
        self.shared_rows_aggregated_calls = pd.merge(self.aggregated_df, self.calls_df.drop(columns=[VAF, ALT_DEPTH]), how="inner",
                                                    on=self.on_list_no_alt, suffixes=("_x", ""))

        # inner join between calls_df and primary_df
        self.shared_rows_primary_calls = pd.merge(self.calls_df.drop(columns=[VAF]), self.primary_df, how="inner", on=self.on_list, suffixes=("_x", ""))


        # rows leftover after we combine shared_rows_aggregated_calls and primary_df and drop duplicates
        self.rows_pulled_from_calls = pd.concat([self.shared_rows_aggregated_calls, self.primary_df]).drop_duplicates(subset=self.on_list_no_alt, keep=False)

        # rows leftover after we combine shared_rows_aggregated_calls and aggregated_df and drop duplicates
        self.imputed_rows = pd.concat([self.aggregated_df, self.shared_rows_aggregated_calls]).drop_duplicates(subset=self.on_list_no_alt, keep=False)

        # number of expected unique loci
        self.n_unique_chr_pos = len(self.test_aggregator.unique_chr_pos) * len(self.test_aggregator.populations)


    def test_df_sizes(self):
        self.assertTrue(len(self.aggregated_df) == self.n_unique_chr_pos,
                        "Size of aggregated_df is not equal to the number of unique chr_pos * number of populations")

        self.assertTrue(len(self.shared_rows_primary_calls) == len(self.primary_df),
                        "All rows in primary_df are not found in calls_df")


        self.assertTrue(len(self.shared_rows_aggregated_primary) == len(self.primary_df),
                        "All rows in primary_df are not in aggregated_df")


    def test_genes_match(self):
        self.assertTrue(all(self.aggregated_df.groupby(CHR_POS)[GENE].nunique().eq(1)) and \
                        (self.aggregated_df.groupby(CHR_POS)[GENE].value_counts().index[0] == self.primary_df.groupby(CHR_POS)[GENE].value_counts().index[0]),
                        "Genes do not match between the aggregated_df and the primary_df for each unique chr_pos")

    def test_imputed_rows(self):
        # number of rows from aggregated_df that have altDepth = 0
        n_zero_altDepth_aggregated = len(self.aggregated_df[self.aggregated_df[ALT_DEPTH] == 0])

        # number of rows that came from primary dataframe that have an altDepth = 0
        n_zero_altDepth_from_primary = len(self.shared_rows_aggregated_primary[self.shared_rows_aggregated_primary[ALT_DEPTH] == 0])

        # the number of rows we should impute across all samples and generate a refDepth for
        n_imputed_zero_altDepth = n_zero_altDepth_aggregated - (len(self.shared_rows_aggregated_calls) - len(self.primary_df))

        self.assertTrue(len(self.imputed_rows) == n_imputed_zero_altDepth,
                        "Incorrect number of imputed rows")


    def test_data_matches(self):

        self.assertTrue(self.shared_rows_aggregated_primary[self.on_list].reset_index(drop=True).equals(self.primary_df[self.on_list].reset_index(drop=True)),
                        "Data does not match between primary_df and aggregated_df")


        self.assertTrue(self.shared_rows_aggregated_primary[self.on_list].reset_index(drop=True).equals(self.primary_df[self.on_list].reset_index(drop=True)),
                        "Data does not match between primary_df and aggregated_df")


if __name__ == '__main__':
    unittest.main()
