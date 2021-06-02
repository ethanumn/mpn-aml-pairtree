from ssm_base_processor import SSM_Base_Processor

class MPN_AML_Processor(SSM_Base_Processor):
    """
    Used to process data collect from AML study into a usable ssm file
    """


    def __init__(self, in_file="", out_file="", wait_to_write_out_file=False):

        super().__init__(in_file, out_file, wait_to_write_out_file)


    def p_names(self):
        """
        Combines the chromosome number stored in the 'seqnames' column (e.g. 'chr1')
        with the starting position for a SNP in that chromosome
        """

        self.out_df[self.COL_NAME] = self.in_df.seqnames + "_" + self.in_df.start.astype(str)


    def p_ids(self):
        """
        Combines the chromosome number stored in the 'seqnames' column (e.g. 'chr1')
        with the starting position for a SNP in that chromosome
        """

        self.out_df[self.COL_ID] = ["s" + str(number) for number in list(range(0, len(self.in_df)))]


    def p_var_reads(self):

        self.out_df[self.COL_VAR_READS] = self.in_df.altDepth


    def p_total_reads(self):

        self.out_df[self.COL_TOTAL_READS] = self.in_df.altDepth + self.in_df.refDepth


    def p_var_read_prob(self):

        self.out_df[self.COL_VAR_READ_PROB] = len(self.in_df) * [0.5] # assumption that everything is diploid
