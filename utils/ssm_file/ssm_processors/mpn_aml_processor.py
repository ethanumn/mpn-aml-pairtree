from ssm_base_processor import SSM_Base_Processor

class MPN_AML_Processor(SSM_Base_Processor):
    """
    Used to process data collect from AML study into a usable ssm file
    """


    def __init__(self, in_file="", out_file="", write_out_file=True, write_out_params=True):

        super().__init__(in_file, out_file, write_out_file, write_out_params)


    def p_names(self):

        self.out_df[self.COL_NAME] = self.in_df[self.CHR_POS]


    def p_ids(self):
        """
        Provides each <chromosome><position> pair with a unique id (r's\d+').
        This column is added to the in_df to aid in writing a params.json file.
        """
        self.in_df[self.COL_ID] = ["s" + str(number) for number in list(range(0, len(self.in_df)))]
        self.out_df[self.COL_ID] = self.in_df[self.COL_ID]


    def p_var_reads(self):

        self.out_df[self.COL_VAR_READS] = self.in_df[self.ALT_DEPTH]


    def p_total_reads(self):

        self.out_df[self.COL_TOTAL_READS] = self.in_df[self.ALT_DEPTH] + self.in_df[self.REF_DEPTH]


    def p_var_read_prob(self):
        
        import re

        self.out_df.loc[self.in_df[self.CHR].str.match('(ch)(.*)((x|y))', flags=re.IGNORECASE), self.COL_VAR_READ_PROB] = 1.0
        self.out_df.loc[self.in_df[self.CHR].str.match('(ch)(.*)(\d+)', flags=re.IGNORECASE), self.COL_VAR_READ_PROB] = 0.5
