import pandas as pd

class MPN_AML_Aggregator:
    """
    A one-off class for aggregating data for the mpn-aml-pairtree analysis
    """
    def __init__(self,
                 master_xls = [],
                 calls_xls = [],
                 populations_xls = [],
                 aggregated_file = [],
                 write_xls_file = True):

        """
        Aims to load in xlsx files, and then kick off preprocessing, processing, and simple verification checks
        """

        self.master_df = self.read_xls_sheet(*master_xls)
        self.calls_df = self.read_xls_sheet(*calls_xls)
        self.populations = self.read_xls_sheet(*populations_xls)[0] # we only want the first column
        self.aggregated_df = pd.DataFrame()
        self.aggregated_file = aggregated_file

        # initialize constants before preprocessing dataframes (or doing anything else for that matter)
        self.init_constants()

        # preprocess dataframes
        self.preprocess_dfs()

        # start processing
        self.process()

        # basic checks to verify aggregation
        self.verify_aggregation()

        # write aggregated and updated dataframe to xls file
        if write_xls_file:

            self.write_xls_sheet(self.aggregated_df, *self.aggregated_file)


    def init_constants(self):
        """
        Initialize constants that are used in aggregating data from different spread sheets
        """
        # column headers used in dataframes

        # for calls df only
        self.SEQNAMES = "seqnames"
        self.START = "start"

        # for master df only
        self.CHR = "Chr"
        self.POSITION = "Position"

        # used by both
        self.CHR_POS = "chr_pos"
        self.REF_DEPTH = "refDepth"
        self.ALT_DEPTH = "altDepth"
        self.SAMPLE_NAMES = "sampleNames"
        self.VAF = "VAF"


        # other constants

        self.SCAN_FILE_EXT = ".RAW.VarScan.txt"



    def read_xls_sheet(self, file_name, sheet_name, header=0):

        return pd.read_excel(file_name, sheet_name, header=header)


    def write_xls_sheet(self, dataframe, file_name, sheet_name):

        if not dataframe.empty:

            dataframe.to_excel(file_name, sheet_name=sheet_name)


    def preprocess_dfs(self):
        """
        Preprocess the different dataframes (e.g. standardize column headers)
        """

        # basic prepocessing to make sure dataframes have same column names, as well as similar formatting of columns
        self.master_df[self.SAMPLE_NAMES] = self.master_df[self.SAMPLE_NAMES].str.replace(self.SCAN_FILE_EXT, "", regex=False)
        self.calls_df[self.SAMPLE_NAMES] = self.calls_df[self.SAMPLE_NAMES].str.replace(self.SCAN_FILE_EXT, "", regex=False)

        # list of columns that each dataframe should have
        self.calls_columns = [self.SEQNAMES, self.START, self.REF_DEPTH, self.ALT_DEPTH, self.SAMPLE_NAMES, self.VAF]
        self.master_columns = [self.CHR, self.POSITION, self.REF_DEPTH, self.ALT_DEPTH, self.SAMPLE_NAMES, self.VAF]
        self.aggregated_columns = [self.CHR, self.POSITION, self.REF_DEPTH, self.ALT_DEPTH, self.SAMPLE_NAMES, self.VAF, self.CHR_POS]

        self.calls_df = self.calls_df[self.calls_columns]
        self.master_df = self.master_df[self.master_columns]

        # rename columns from calls dataframe
        self.calls_df = self.calls_df.rename(columns={self.SEQNAMES: self.CHR, self.START: self.POSITION})

        # add <chromosome><position> column which will help us compile a proper aggregated dataframe
        self.master_df[self.CHR_POS] = self.master_df[self.CHR] + "_" + self.master_df[self.POSITION].apply(str)
        self.calls_df[self.CHR_POS] = self.calls_df[self.CHR] + "_" + self.calls_df[self.POSITION].apply(str)

        self.unique_chr_pos = self.master_df[self.CHR_POS].unique() # obtain all unique chromosome + position pairs


    def init_aggregated_df(self):
        """
        Create an empty dataframe containing a unique <chromosome><position> found in the master
        dataframe for each sample (n_samples * len(self.unique_chr_pos))
        """

        for pop in self.populations:

            self.aggregated_df = self.aggregated_df.append(pd.DataFrame({

                self.CHR          : pd.Series([chr_pos.split("_")[0] for chr_pos in self.unique_chr_pos], dtype="object"),
                self.POSITION     : pd.Series([chr_pos.split("_")[1] for chr_pos in self.unique_chr_pos], dtype="int64"),
                self.CHR_POS      : pd.Series(self.unique_chr_pos, dtype="object"),
                self.REF_DEPTH    : pd.Series([], dtype="int64"),
                self.ALT_DEPTH    : pd.Series([], dtype="int64"),
                self.SAMPLE_NAMES : pd.Series([pop]*len(self.unique_chr_pos), dtype="object"),
                self.VAF          : pd.Series([], dtype="float64"),

            }))


    def impute_missing_values(self):
        import math

        # impute ref depth values for
        for chr_pos in self.unique_chr_pos:

            # obtain every samples <chromosome><position> pair that does not have a NaN refDepth value
            chr_pos_df = self.aggregated_df.loc[(self.aggregated_df[self.CHR_POS] == chr_pos)
                                           & self.aggregated_df[self.REF_DEPTH].notnull()]

            # overwrite NaN values with the average total reads across all samples for a <chromosome><position> pair
            self.aggregated_df.loc[(self.aggregated_df[self.CHR_POS] == chr_pos)
                                   & self.aggregated_df[self.REF_DEPTH].isnull(), self.REF_DEPTH] = \
                math.floor((chr_pos_df[self.REF_DEPTH].sum() + chr_pos_df[self.ALT_DEPTH].sum()) / len(chr_pos_df))


    def process(self):
        """
        Initializes dataframe for aggegration, then merges all dataframes
        """

        on_list = [self.CHR_POS, self.POSITION, self.CHR, self.SAMPLE_NAMES] # columns to use for join

        self.init_aggregated_df()

        # join the master dataframe with the aggregated dataframe
        self.aggregated_df = self.aggregated_df.merge(self.master_df, how="left", on=on_list, suffixes=("_x", ""))

        # join the calls dataframe with the aggregated dataframe
        self.aggregated_df = self.aggregated_df.merge(self.calls_df, how="left", on=on_list, suffixes=("_x", ""))

        # fill all NaN values in altDepth column as 0
        self.aggregated_df[self.ALT_DEPTH] = self.aggregated_df[self.ALT_DEPTH].fillna(0)

        # fill all NaN values in VAF column as 0
        self.aggregated_df[self.VAF] = self.aggregated_df[self.VAF].fillna(0)


        # impute ref depth values for missing variants
        self.impute_missing_values()

        # get rid of any extra columns from the merges
        self.aggregated_df = self.aggregated_df[self.aggregated_columns]

        # need to reset column type since joins will change the underlying type (and so we"ll fail our checks)
        self.aggregated_df[self.POSITION] = self.aggregated_df[self.POSITION].astype("int64")
        self.aggregated_df[self.REF_DEPTH] = self.aggregated_df[self.REF_DEPTH].astype("int64")
        self.aggregated_df[self.ALT_DEPTH] = self.aggregated_df[self.ALT_DEPTH].astype("int64")


    def verify_aggregation(self):


        # list of columns to join on - do not use VAF or any floating point value since the changes in precision
        # do not allow for a proper join
        on_list = [self.CHR_POS, self.CHR, self.POSITION, self.SAMPLE_NAMES, self.ALT_DEPTH, self.REF_DEPTH]


        # inner join combinations
        shared_rows_aggregate_master = pd.merge(self.master_df, self.aggregated_df, how="inner", on=on_list)
        shared_rows_aggregate_calls = pd.merge(self.calls_df, self.aggregated_df, how="inner", on=on_list)
        shared_rows_master_calls = pd.merge(self.master_df, self.calls_df, how="inner", on=on_list)

        # number of possible imputed rows
        zero_altDepth_aggregate = len(self.aggregated_df[self.aggregated_df[self.ALT_DEPTH] == 0])

        # number of rows that came from calls dataframe that have an altDepth = 0
        zero_altDepth_from_calls = len(shared_rows_aggregate_calls[shared_rows_aggregate_calls[self.ALT_DEPTH] == 0])

        # number of rows that came from master dataframe that have an altDepth = 0
        zero_altDepth_from_master = len(shared_rows_aggregate_master[shared_rows_aggregate_master[self.ALT_DEPTH] == 0])

        n_unique_chr_pos = len(self.unique_chr_pos) * len(self.populations)

        n_imputed_zero_altDepth = zero_altDepth_aggregate - (zero_altDepth_from_calls + zero_altDepth_from_master)

        # the number of values we should impute are the total number of <chromosome><position> pairs
        # minus the number of row we found in the master df, and calls df, while accounting for duplicates
        no_entry_master_or_calls = n_unique_chr_pos - (len(shared_rows_aggregate_master) +
                                                      len(shared_rows_aggregate_calls) -
                                                       len(shared_rows_master_calls))

        print("number of rows in aggregate / ((# of samples) * (# of unique <chromosome><position> pairs)) = %d/%d"
             % (len(self.aggregated_df), len(self.master_df[self.CHR_POS].unique()) * len(self.populations)))

        print("number of matching rows between aggregate and master / number of rows in master = %d/%d"
              % (len(shared_rows_aggregate_master), len(self.master_df)))

        print("number of matching rows between aggregate and calls = %d"
              % (len(shared_rows_aggregate_calls)))

        print("number of <chromosome><position> pairs with altDepth = 0 / remaining unmatched rows (from either calls or master) = %d/%d"
             % (n_imputed_zero_altDepth, no_entry_master_or_calls))
