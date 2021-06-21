import sys, os

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from mpn_aml_columns import *
from ssm_columns import *
from ssm_base_processor import SSM_Base_Processor

class MPN_AML_Processor(SSM_Base_Processor):
    """
    Used to process data collect from AML study into a usable ssm file
    """


    def __init__(self, in_file="", out_file="", write_out_file=True, write_out_params=True):

        super().__init__(in_file, out_file, write_out_file, write_out_params)


    def format_out_df(self):
        """
        Override to translate the processed_df into a dataframe
        which will be written out as a .ssm file
        """
        import pandas as pd

        no_brackets_array_string = lambda arr: [", ".join(map(str, entry)) for entry in arr]

        # initialize out_df
        self.out_df = pd.DataFrame()

        # set name column <gene>_<position> which we assume to be unique for a given <chromosome><position>
        self.out_df[COL_NAME] = self.processed_df.groupby(COL_NAME)[COL_NAME].agg(pd.Series.mode).values

        # set var_reads
        self.out_df[COL_VAR_READS] = no_brackets_array_string(self.processed_df.groupby(COL_NAME)[COL_VAR_READS].apply(list).values)

        # set total_reads
        self.out_df[COL_TOTAL_READS] = no_brackets_array_string(self.processed_df.groupby(COL_NAME)[COL_TOTAL_READS].apply(list).values)

        # set var_read_prob
        self.out_df[COL_VAR_READ_PROB] = no_brackets_array_string(self.processed_df.groupby(COL_NAME)[COL_VAR_READ_PROB].apply(list).values)

        # set name
        # Provides each <chromosome><position> pair with a unique id (r's\d+').
        # We're purposely doing this after all of the other columns have been created.
        self.out_df[COL_ID] = ["s" + str(number) for number in list(range(0, len(self.out_df)))]


    def p_names(self):

        self.processed_df[COL_NAME] = self.in_df[GENE] + "_" + self.in_df[POSITION].apply(str)


    def p_var_reads(self):

        self.processed_df[COL_VAR_READS] = self.in_df[ALT_DEPTH]


    def p_total_reads(self):

        self.processed_df[COL_TOTAL_READS] = self.in_df[ALT_DEPTH] + self.in_df[REF_DEPTH]


    def p_var_read_prob(self):

        import re

        self.processed_df.loc[self.in_df[CHR].str.match('(ch)(.*)((x|y))', flags=re.IGNORECASE), COL_VAR_READ_PROB] = 1.0
        self.processed_df.loc[self.in_df[CHR].str.match('(ch)(.*)(\d+)', flags=re.IGNORECASE), COL_VAR_READ_PROB] = 0.5
