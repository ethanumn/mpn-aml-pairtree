import unittest
import numpy as np
import pandas as pd
import os, sys

sys.path.append("../ssm_processors")

from mpn_aml_processor import MPN_AML_Processor

class MPN_AML_Processor_Tests(unittest.TestCase):

    def setUp(self):
        self.in_file = os.environ["DATA_DIR"] + "/" + "MATS08_sq.xlsx"
        self.out_file = os.environ["DATA_DIR"] + "/" + "test.ssm"
        self.wait_to_write_out_file = True

        self.ssm_obj = MPN_AML_Processor(self.in_file, self.out_file, self.wait_to_write_out_file)


    def test_size(self):
        self.assertEqual(len(self.ssm_obj.out_df), len(self.ssm_obj.in_df),
                        'Incorrect output size')


    def test_name_convention(self):
        self.assertTrue(all(self.ssm_obj.out_df.id.str.match(r's\d+') == True),
                         'Incorrect id naming convention')


    def test_vaf(self):
        self.assertTrue(
            all(
                np.isclose(
                    self.ssm_obj.in_df.VAF,
                    (self.ssm_obj.out_df.var_reads / self.ssm_obj.out_df.total_reads),
                    atol=1e-10
                )
                == True
            ),
            'Incorrect variant allele frequency calculation'
        )


    def test_var_reads(self):
        self.assertTrue(pd.Series.equals(self.ssm_obj.in_df.altDepth, self.ssm_obj.out_df.var_reads),
                       'Incorrect variant read')


    def test_total_reads(self):
        self.assertTrue(pd.Series.equals(self.ssm_obj.in_df.altDepth+self.ssm_obj.in_df.refDepth,
                                         self.ssm_obj.out_df.total_reads),
                       'Incorrect total reads calculation')


if __name__ == '__main__':
    unittest.main()
