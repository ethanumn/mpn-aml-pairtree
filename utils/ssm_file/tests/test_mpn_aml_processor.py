import unittest
import os, sys
import argparse

import numpy as np
import pandas as pd

sys.path.append(os.environ["UTILS_DIR"] + "/common")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ssm_processors'))

from mpn_aml_columns import *
from ssm_columns import *
from mpn_aml_processor import MPN_AML_Processor

class MPN_AML_Processor_Tests(unittest.TestCase):
    """
    Test cases for MPN_AML_Processor for converting custom file format to .ssm.
    use 'python3 test_mpn_aml_processor.py' to run the test suite
    """
    def setUp(self):

        # set up using example in-file
        self.in_file = os.environ["DATA_DIR"] + "/example/results/" + "example.aggregated.xlsx" # change this to test with a different file
        #self.in_file = os.environ["DATA_DIR"] + "/mpn_aml_test/results/" + "test.aggregated.xlsx" # change this to test with a different file

        write_out_file = False
        write_out_params = False

        self.test_processor = MPN_AML_Processor(self.in_file, "", write_out_file, write_out_params)

        # rename for simplification of statements
        self.processed_df = self.test_processor.processed_df
        self.in_df = self.test_processor.in_df
        self.out_df = self.test_processor.out_df

    def test_processed_size(self):
        # in_df is aggregated file - so the size should be the same as processed_df but the contents will differ
        self.assertEqual(len(self.processed_df), len(self.in_df),
                        'Processed dataframe has incorrect output size')


    def test_processed_vaf(self):
        # the VAF of the values in the in_df (that existed the primary xls (that were not imputed),
        # should be the same (close) to the calculation we get using the processed_df
        self.assertTrue(
            all(
                np.isclose(
                    self.in_df.loc[self.in_df[ALT_DEPTH] != 0, VAF],
                    (self.processed_df.loc[self.processed_df[COL_VAR_READS] != 0, COL_VAR_READS] /
                     self.processed_df.loc[self.processed_df[COL_VAR_READS] != 0, COL_TOTAL_READS]),
                    atol=1e-10
                )
                == True
            ),
            'Processed dataframe has incorrect variant allele frequency calculation'
        )


    def test_processed_var_reads(self):
        # the altDepth from in_df should be the same as the var_reads in processed_df
        self.assertTrue(pd.Series.equals(self.in_df[ALT_DEPTH], self.processed_df[COL_VAR_READS]),
                       'Processed dataframe has incorrect variant read')


    def test_processed_total_reads(self):
        # total_reads in processed_df should be equal to refDepth+altDepth from in_df
        self.assertTrue(pd.Series.equals(self.in_df[ALT_DEPTH]+self.in_df[REF_DEPTH],
                                         self.processed_df[COL_TOTAL_READS]),
                       'Processed dataframe has incorrect total reads calculation')


    def test_out_id_convention(self):
        # id convention should match regular expression
        self.assertTrue(all(self.out_df[COL_ID].str.match(r's\d+') == True),
                         'Processed dataframe has incorrect id convention')


    def test_out_name(self):
        # there should be one entry for each of the unique name's in the processed_df in the out_df
        self.assertTrue(sorted(self.processed_df[COL_NAME].unique()) == sorted(self.out_df[COL_NAME].values),
                       'Out dataframe has incorrect values in the name column')


    def test_out_var_reads(self):
        # the var_read groups in the out_df should be the same as grouping all of the altDepth reads per unique <chromosome><position> (when sorted that is)
        self.assertTrue(sorted(np.array(self.in_df.groupby(CHR_POS)[ALT_DEPTH].apply(list).values).transpose()) == sorted([[int(i) for i in arr.split(",")] for arr in self.out_df[COL_VAR_READS].values]),
                       'Processed dataframe has incorrect variant read')


    def test_out_vector_sizes(self):
        # there should be the same amount of (var_reads, total_reads, var_read_prob)
        # in each grouping of <chromosome><position> in the out_df as there are samples in the in_df
        self.assertTrue(all(len(self.in_df[SAMPLE_NAMES].unique()) == self.out_df[COL_VAR_READS].apply(lambda arr: len(arr.split(","))).values),
                       'The size of the vectors in var_reads fields is incorrect')
        self.assertTrue(all(len(self.in_df[SAMPLE_NAMES].unique()) == self.out_df[COL_TOTAL_READS].apply(lambda arr: len(arr.split(","))).values),
                      'The size of the vectors in total_reads fields is incorrect')
        self.assertTrue(all(len(self.in_df[SAMPLE_NAMES].unique()) == self.out_df[COL_VAR_READ_PROB].apply(lambda arr: len(arr.split(","))).values),
                     'The size of the vectors in var_read_prob fields is incorrect')




if __name__ == '__main__':
    unittest.main()
