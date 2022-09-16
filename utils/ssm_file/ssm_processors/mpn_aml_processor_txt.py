import sys, os
import numpy as np
import pandas as pd

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from mpn_aml_columns_txt import *
from ssm_columns import *
from ssm_base_processor import SSM_Base_Processor

class MPN_AML_Processor_Txt(SSM_Base_Processor):
    """
    Used to process data collect from AML study into a usable ssm file
    """


    def __init__(self, in_file="", out_file="", write_out_file=True, write_out_params=True):

        super().__init__(in_file, out_file, write_out_file, write_out_params, SAMPLEA, sort_samples=True)
    

    def format_out_df(self):
        """
        Override to translate the processed_df into a dataframe
        which will be written out as a .ssm file
        """

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

        # make variants be listed in ascending order (by chromosome and position)
        idx1 = self.out_df[COL_NAME].apply(
                lambda n: 
                    int(
                        n.lstrip("chr").split("_")[0].replace("X", "98").replace("Y","99") 
                    )
            )
            
        idx2 = self.out_df[COL_NAME].apply(
                lambda n: 
                    int(
                        n.split("_")[1])
                    )
                    
        self.out_df = self.out_df.set_index([idx1, idx2])
        self.out_df = self.out_df.sort_index()
        # set name
        # Provides each <chromosome><position> pair with a unique id (r's\d+').
        # We're purposely doing this after all of the other columns have been created.
        self.out_df[COL_ID] = ["s" + str(number) for number in list(range(0, len(self.out_df)))]
        

    def p_df_sort(self):
        # sort dataframe by sample 
        self.in_df = self.in_df.sort_values(by=[SAMPLEA])
        
        
    def p_names(self):

        self.processed_df[COL_NAME] = self.in_df[CHR].astype(str) + "_" + self.in_df[START].astype(str)


    def p_var_reads(self):

        self.processed_df[COL_VAR_READS] = self.in_df[VAR_READS]
        self.processed_df[COL_VAR_READS] = self.processed_df[COL_VAR_READS].fillna(0).astype(int)


    def p_total_reads(self):

        self.processed_df[COL_TOTAL_READS] = self.in_df[TOTAL_READS]
        self.processed_df[COL_TOTAL_READS] = self.processed_df[COL_TOTAL_READS].fillna(1).astype(int)
        self.processed_df.loc[self.processed_df[COL_TOTAL_READS] < 1, COL_TOTAL_READS] = 1
        
        
    def p_var_read_prob(self):

        import re
        
        # if it's a sex chromosome, then it should have a var_read_prob = 1
        self.processed_df.loc[self.in_df[CHR].str.match('(chr)(.*)((x|y))', flags=re.IGNORECASE), COL_VAR_READ_PROB] = 1.0
        
        # if it's an autosome, it should it should have a var_read_prob of M/N
        self.processed_df.loc[self.in_df[CHR].str.match('(chr)(.*)(\d+)', flags=re.IGNORECASE), COL_VAR_READ_PROB] = \
                round(1 / (2 + (self.in_df[COPY_NUMBER].astype(float) - 2)), 3)
                
        # if we have a copy number of 1, we're assuming a LOH of the reference allele
        self.processed_df.loc[self.in_df[COPY_NUMBER] == 1, COL_VAR_READ_PROB] = 1.0
        
        # if we have a copy number of 0, we'll assume it's not present
        self.processed_df.loc[self.in_df[COPY_NUMBER] == 0, COL_VAR_READ_PROB] = np.nan
        
        print("# rows BEFORE dropping nans: %d" % len(self.processed_df))
        
        # we can't use variants with a var_read_prob of 0 or nan
        self.processed_df = self.processed_df.dropna()
        
        print("# rows AFTER dropping nans: %d" % len(self.processed_df))
