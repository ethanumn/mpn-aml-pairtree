import pandas as pd
import sys, os


sys.path.append(os.environ["UTILS_DIR"] + "/pdf_templates")
sys.path.append(os.environ["UTILS_DIR"] + "/common")

from mpn_aml_columns import *
from utils.verify_aggregation import verify_aggregation

# impute techniques
IMPUTE_AVG = "AVG"
IMPUTE_ZERO = "ZERO"


class MPN_AML_Aggregator:
    """
    A one-off class for aggregating data for the mpn-aml-pairtree analysis
    """
    def __init__(self,
                 primary_xls = [],
                 calls_xls = [],
                 populations_xls = [],
                 aggregated_xls = [],
                 metrics_file = "",
                 write_xls_file = True,
                 impute_technique=IMPUTE_ZERO):

        """
        Aims to load in xlsx files, and then kick off preprocessing, processing, and simple verification checks
        """

        self.primary_df = self.read_xls_sheet(*primary_xls)
        self.calls_df = self.read_xls_sheet(*calls_xls)
        self.populations = self.read_xls_sheet(*populations_xls)[0] # we only want the first column
        self.aggregated_df = pd.DataFrame()
        self.aggregated_xls = aggregated_xls
        self.metrics_file = metrics_file
        self.impute_technique = impute_technique

        # initialize constants before preprocessing dataframes (or doing anything else for that matter)
        self.init_constants()

        # preprocess dataframes
        self.preprocess_dfs()

        # start processing
        self.process()

        # basic checks to verify aggregation
        verify_aggregation(self.metrics_file,
                           self.aggregated_df,
                           self.primary_df,
                           self.calls_df,
                           self.populations,
                           self.unique_chr_pos,
                           self.aggregated_columns,
                           primary_xls,
                           calls_xls,
                           populations_xls,
                           aggregated_xls)

        # write aggregated and updated dataframe to xls file
        if write_xls_file:

            self.write_xls_sheet(self.aggregated_df, *self.aggregated_xls)


    def init_constants(self):
        """
        Initialize constants that are used in aggregating data from different spread sheets
        """

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
        self.primary_df[SAMPLE_NAMES] = self.primary_df[SAMPLE_NAMES].str.replace(self.SCAN_FILE_EXT, "", regex=False)
        self.calls_df[SAMPLE_NAMES] = self.calls_df[SAMPLE_NAMES].str.replace(self.SCAN_FILE_EXT, "", regex=False)

        # list of columns that each dataframe should have
        self.calls_columns = [SEQNAMES, START, REF_DEPTH, ALT_DEPTH, SAMPLE_NAMES, VAF]
        self.primary_columns = [CHR, POSITION, REF_DEPTH, ALT_DEPTH, SAMPLE_NAMES, VAF, GENE]
        self.aggregated_columns = [CHR, POSITION, REF_DEPTH, ALT_DEPTH, SAMPLE_NAMES, VAF, GENE, CHR_POS]

        self.calls_df = self.calls_df[self.calls_columns]
        self.primary_df = self.primary_df[self.primary_columns]

        # rename columns from calls dataframe
        self.calls_df = self.calls_df.rename(columns={SEQNAMES: CHR, START: POSITION})

        # add <chromosome><position> column which will help us compile a proper aggregated dataframe
        self.primary_df[CHR_POS] = self.primary_df[CHR] + "_" + self.primary_df[POSITION].apply(str)
        self.calls_df[CHR_POS] = self.calls_df[CHR] + "_" + self.calls_df[POSITION].apply(str)

        # sort the dataframe by chromosome number by creating CHR_NUM column which contains chromsome number, then drop that column after sorting
        self.primary_df[CHR_NUM] = self.primary_df[CHR].str.extract("(\d+)", expand=False).astype(int, errors = "ignore").fillna(0).astype(int)
        self.primary_df = self.primary_df.sort_values(by=[CHR_NUM, POSITION]).drop(columns=[CHR_NUM])

        self.calls_df[CHR_NUM] = self.calls_df[CHR].str.extract("(\d+)", expand=False).astype("int64", errors = "ignore").fillna(0).astype(int)
        self.calls_df = self.calls_df.sort_values(by=[CHR_NUM, POSITION]).drop(columns=[CHR_NUM])

        self.unique_chr_pos = self.primary_df[CHR_POS].unique() # obtain all unique chromosome + position pairs


    def init_aggregated_df(self):
        """
        Create an empty dataframe containing a unique <chromosome><position> found in the primary
        dataframe for each sample (n_samples * len(self.unique_chr_pos))
        """

        self.aggregated_df = pd.DataFrame()

        for pop in self.populations:

            self.aggregated_df = self.aggregated_df.append(pd.DataFrame({

                CHR          : pd.Series([chr_pos.split("_")[0] for chr_pos in self.unique_chr_pos], dtype="object"),
                POSITION     : pd.Series([chr_pos.split("_")[1] for chr_pos in self.unique_chr_pos], dtype="int64"),
                CHR_POS      : pd.Series(self.unique_chr_pos, dtype="object"),
                REF_DEPTH    : pd.Series([], dtype="int64"),
                ALT_DEPTH    : pd.Series([], dtype="int64"),
                SAMPLE_NAMES : pd.Series([pop]*len(self.unique_chr_pos), dtype="object"),
                GENE         : pd.Series([], dtype="object"),
                VAF          : pd.Series([], dtype="float64")

            }))

        # sort the dataframe by chromosome number by creating CHR_NUM column which contains chromsome number, then drop that column after sorting
        self.aggregated_df[CHR_NUM] = self.aggregated_df[CHR].str.extract("(\d+)", expand=False).astype(int, errors = "ignore").fillna(0).astype(int)
        self.aggregated_df = self.aggregated_df.sort_values(by=[CHR_NUM, POSITION]).drop(columns=[CHR_NUM])


    def impute_missing_values(self):

        if self.impute_technique == IMPUTE_AVG:

            import math

            # impute ref depth values for
            for chr_pos in self.unique_chr_pos:

                # obtain every samples <chromosome><position> pair that does not have a NaN refDepth value
                chr_pos_df = self.primary_df.loc[(self.primary_df[CHR_POS] == chr_pos)
                                               & self.primary_df[REF_DEPTH].notnull()]

                # overwrite NaN values with the average total reads across all samples for a <chromosome><position> pair
                self.aggregated_df.loc[(self.aggregated_df[CHR_POS] == chr_pos)
                                       & self.aggregated_df[REF_DEPTH].isnull(), REF_DEPTH] = \
                    math.floor((chr_pos_df[REF_DEPTH].sum() + chr_pos_df[ALT_DEPTH].sum()) / len(chr_pos_df))

        elif self.impute_technique == IMPUTE_ZERO:

            self.aggregated_df[REF_DEPTH] = self.aggregated_df[REF_DEPTH].fillna(1)


    def process(self):
        """
        Initializes dataframe for aggegration, then merges all dataframes
        """

        on_list = [CHR_POS, POSITION, CHR, SAMPLE_NAMES] # columns to use for join

        self.init_aggregated_df()

        # join the primary dataframe with the aggregated dataframe
        self.aggregated_df = self.aggregated_df.merge(self.primary_df, how="left", on=on_list, suffixes=("_x", ""))

        # join the calls dataframe with the aggregated dataframe
        self.aggregated_df = self.aggregated_df.merge(self.calls_df.drop(columns=[ALT_DEPTH, VAF]), how="left", on=on_list, suffixes=("_x", ""))

        # fill all NaN values in altDepth column as 0
        # we're doing this here because we want to set variant reads to 0 for all <chromosome><position> pairs for a sample that weren't in the primary spreadsheet
        self.aggregated_df[ALT_DEPTH] = self.aggregated_df[ALT_DEPTH].fillna(0)

        # fill all NaN values in VAF column as 0
        self.aggregated_df[VAF] = self.aggregated_df[VAF].fillna(0)

        # fill gene column with the most common gene that the matching <chromosome><position> pairs have
        self.aggregated_df[GENE] =  self.aggregated_df.groupby(CHR_POS)[GENE].transform(lambda grp: grp.fillna(grp[grp.notnull()].mode()[0]))

        # impute ref depth values for missing variants
        self.impute_missing_values()

        # get rid of any extra columns from the merges
        self.aggregated_df = self.aggregated_df[self.aggregated_columns]

        # need to reset column type since joins will change the underlying type (and so we"ll fail our checks)
        self.aggregated_df[POSITION] = self.aggregated_df[POSITION].astype("int64")
        self.aggregated_df[REF_DEPTH] = self.aggregated_df[REF_DEPTH].astype("int64")
        self.aggregated_df[ALT_DEPTH] = self.aggregated_df[ALT_DEPTH].astype("int64")
