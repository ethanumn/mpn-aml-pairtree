import unittest
import numpy as np
import pandas as pd
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ssm_processors'))

from mpn_aml_processor import MPN_AML_Processor

class MPN_AML_Processor_Tests(unittest.TestCase):
    """
    Test cases for MPN_AML_Processor for converting custom file format to .ssm.
    use 'python3 test_mpn_aml_processor.py' to run the test suite
    """
    def setUp(self):

        # set up using example in-file
        self.in_file = os.environ["DATA_DIR"] + "/example/" + "example.aggregated.xlsx"

        write_out_file = False
        write_out_params = False

        self.test_processor = MPN_AML_Processor(self.in_file, "", write_out_file, write_out_params)


    def test_size(self):
        self.assertEqual(len(self.test_processor.out_df), len(self.test_processor.in_df),
                        'Incorrect output size')


    def test_name_convention(self):
        self.assertTrue(all(self.test_processor.out_df.id.str.match(r's\d+') == True),
                         'Incorrect id naming convention')


    def test_vaf(self):
        self.assertTrue(
            all(
                np.isclose(
                    self.test_processor.in_df.VAF,
                    (self.test_processor.out_df.var_reads / self.test_processor.out_df.total_reads),
                    atol=1e-10
                )
                == True
            ),
            'Incorrect variant allele frequency calculation'
        )


    def test_var_reads(self):
        self.assertTrue(pd.Series.equals(self.test_processor.in_df.altDepth, self.test_processor.out_df.var_reads),
                       'Incorrect variant read')


    def test_total_reads(self):
        self.assertTrue(pd.Series.equals(self.test_processor.in_df.altDepth+self.test_processor.in_df.refDepth,
                                         self.test_processor.out_df.total_reads),
                       'Incorrect total reads calculation')


if __name__ == '__main__':
    unittest.main()
