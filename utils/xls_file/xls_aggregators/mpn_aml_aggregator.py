import pandas as pd
import sys, os

sys.path.append(os.environ["UTILS_DIR"] + "/pdf_templates")

from mpn_aml_metrics_pdf import MPN_AML_METRICS_PDF

class MPN_AML_Aggregator:
    """
    A one-off class for aggregating data for the mpn-aml-pairtree analysis
    """
    def __init__(self,
                 master_xls = [],
                 calls_xls = [],
                 populations_xls = [],
                 aggregated_xls = [],
                 metrics_file = "",
                 write_xls_file = True):

        """
        Aims to load in xlsx files, and then kick off preprocessing, processing, and simple verification checks
        """

        self.master_df = self.read_xls_sheet(*master_xls)
        self.calls_df = self.read_xls_sheet(*calls_xls)
        self.populations = self.read_xls_sheet(*populations_xls)[0] # we only want the first column
        self.aggregated_df = pd.DataFrame()
        self.aggregated_xls = aggregated_xls
        self.metrics_file = metrics_file

        # initialize constants before preprocessing dataframes (or doing anything else for that matter)
        self.init_constants()

        # preprocess dataframes
        self.preprocess_dfs()

        # start processing
        self.process()

        # basic checks to verify aggregation
        self.verify_aggregation(master_xls, calls_xls, populations_xls, aggregated_xls)

        # write aggregated and updated dataframe to xls file
        if write_xls_file:

            self.write_xls_sheet(self.aggregated_df, *self.aggregated_xls)


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
        self.GENE = "gene"


        # other constants
        self.CHR_NUM = "chr_num"

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
        self.master_columns = [self.CHR, self.POSITION, self.REF_DEPTH, self.ALT_DEPTH, self.SAMPLE_NAMES, self.VAF, self.GENE]
        self.aggregated_columns = [self.CHR, self.POSITION, self.REF_DEPTH, self.ALT_DEPTH, self.SAMPLE_NAMES, self.VAF, self.GENE, self.CHR_POS]

        self.calls_df = self.calls_df[self.calls_columns]
        self.master_df = self.master_df[self.master_columns]

        # rename columns from calls dataframe
        self.calls_df = self.calls_df.rename(columns={self.SEQNAMES: self.CHR, self.START: self.POSITION})

        # add <chromosome><position> column which will help us compile a proper aggregated dataframe
        self.master_df[self.CHR_POS] = self.master_df[self.CHR] + "_" + self.master_df[self.POSITION].apply(str)
        self.calls_df[self.CHR_POS] = self.calls_df[self.CHR] + "_" + self.calls_df[self.POSITION].apply(str)

        # sort the dataframe by chromosome number by creating CHR_NUM column which contains chromsome number, then drop that column after sorting
        self.master_df[self.CHR_NUM] = self.master_df[self.CHR].str.extract("(\d+)", expand=False).astype(int, errors = "ignore").fillna(0).astype(int)
        self.master_df = self.master_df.sort_values(by=[self.CHR_NUM, self.POSITION]).drop(columns=[self.CHR_NUM])

        self.calls_df[self.CHR_NUM] = self.calls_df[self.CHR].str.extract("(\d+)", expand=False).astype("int64", errors = "ignore").fillna(0).astype(int)
        self.calls_df = self.calls_df.sort_values(by=[self.CHR_NUM, self.POSITION]).drop(columns=[self.CHR_NUM])

        self.unique_chr_pos = self.master_df[self.CHR_POS].unique() # obtain all unique chromosome + position pairs


    def init_aggregated_df(self):
        """
        Create an empty dataframe containing a unique <chromosome><position> found in the master
        dataframe for each sample (n_samples * len(self.unique_chr_pos))
        """

        self.aggregated_df = pd.DataFrame()

        for pop in self.populations:

            self.aggregated_df = self.aggregated_df.append(pd.DataFrame({

                self.CHR          : pd.Series([chr_pos.split("_")[0] for chr_pos in self.unique_chr_pos], dtype="object"),
                self.POSITION     : pd.Series([chr_pos.split("_")[1] for chr_pos in self.unique_chr_pos], dtype="int64"),
                self.CHR_POS      : pd.Series(self.unique_chr_pos, dtype="object"),
                self.REF_DEPTH    : pd.Series([], dtype="int64"),
                self.ALT_DEPTH    : pd.Series([], dtype="int64"),
                self.SAMPLE_NAMES : pd.Series([pop]*len(self.unique_chr_pos), dtype="object"),
                self.GENE         : pd.Series([], dtype="object"),
                self.VAF          : pd.Series([], dtype="float64")

            }))

        # sort the dataframe by chromosome number by creating CHR_NUM column which contains chromsome number, then drop that column after sorting
        self.aggregated_df[self.CHR_NUM] = self.aggregated_df[self.CHR].str.extract("(\d+)", expand=False).astype(int, errors = "ignore").fillna(0).astype(int)
        self.aggregated_df = self.aggregated_df.sort_values(by=[self.CHR_NUM, self.POSITION]).drop(columns=[self.CHR_NUM])


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
        self.aggregated_df = self.aggregated_df.merge(self.calls_df.drop(columns=[self.ALT_DEPTH]), how="left", on=on_list, suffixes=("_x", ""))

        # fill all NaN values in altDepth column as 0
        # we're doing this here because we want to set variant reads to 0 for all <chromosome><position> pairs for a sample that weren't in the master spreadsheet
        self.aggregated_df[self.ALT_DEPTH] = self.aggregated_df[self.ALT_DEPTH].fillna(0)

        # fill all NaN values in VAF column as 0
        self.aggregated_df[self.VAF] = self.aggregated_df[self.VAF].fillna(0)

        # fill gene column with the most common gene that the matching <chromosome><position> pairs have
        self.aggregated_df[self.GENE] =  self.aggregated_df.groupby(self.CHR_POS)[self.GENE].transform(lambda grp: grp.fillna(grp[grp.notnull()].mode()[0]))

        # impute ref depth values for missing variants
        self.impute_missing_values()

        # get rid of any extra columns from the merges
        self.aggregated_df = self.aggregated_df[self.aggregated_columns]

        # need to reset column type since joins will change the underlying type (and so we"ll fail our checks)
        self.aggregated_df[self.POSITION] = self.aggregated_df[self.POSITION].astype("int64")
        self.aggregated_df[self.REF_DEPTH] = self.aggregated_df[self.REF_DEPTH].astype("int64")
        self.aggregated_df[self.ALT_DEPTH] = self.aggregated_df[self.ALT_DEPTH].astype("int64")


    def verify_aggregation(self, master_xls="", calls_xls="", populations_xls="", aggregated_xls=""):
        """
        Verify aggregation by producing a pdf that contains the following:
            - Verification of the number of rows in the different dataframes and some unique corner cases
            - Plot/Table of imputed refDepth values for any unique <chromosome><position>
              that should have been in the master xls but was not
            - Plot of VAF of each <chromosome><position> pair per sample
            - Plot/Table of each <chromosome><position> pair pulled from the calls xls per sample
        """
        from tqdm import tqdm


        # list of columns to join on - do not use VAF or any floating point value since the changes in precision
        # do not allow for a proper join. DO NOT use altDepth when comparing between calls_df and aggregate_df because we've zero'd out
        # the altDepth of those rows we pulled from calls_df into the aggregate_df
        on_list = [self.CHR_POS, self.CHR, self.POSITION, self.SAMPLE_NAMES, self.ALT_DEPTH, self.REF_DEPTH]
        on_list_no_alt = [self.CHR_POS, self.CHR, self.POSITION, self.SAMPLE_NAMES, self.REF_DEPTH]
        table_columns = [self.SAMPLE_NAMES, self.CHR_POS, self.POSITION, self.ALT_DEPTH, self.REF_DEPTH]

        # inner join between master_df and aggregate_df
        shared_rows_aggregate_master = pd.merge(self.master_df, self.aggregated_df, how="inner", on=on_list, suffixes=("_x", ""))

        # inner join between calls_df and master_df
        shared_rows_master_calls = pd.merge(self.calls_df.drop(columns=[self.VAF]), self.master_df, how="inner", on=on_list, suffixes=("_x", ""))

        # inner join between aggregate_df and calls_df
        shared_rows_aggregate_calls = pd.merge(self.aggregated_df, self.calls_df.drop(columns=[self.VAF]), how="inner",
                                                    on=on_list_no_alt, suffixes=("_x", "")).dropna()

        # rows that have ref depth imputed
        imputed_rows = pd.concat([self.aggregated_df, shared_rows_aggregate_calls]).drop_duplicates(subset=on_list_no_alt, keep=False)

        # number of possible imputed rows
        n_zero_altDepth_aggregate = len(self.aggregated_df[self.aggregated_df[self.ALT_DEPTH] == 0])

        # number of rows that came from calls dataframe that have an altDepth = 0
        n_zero_altDepth_from_calls = len(shared_rows_aggregate_calls[shared_rows_aggregate_calls[self.ALT_DEPTH] == 0])

        # number of rows that came from master dataframe that have an altDepth = 0
        n_zero_altDepth_from_master = len(shared_rows_aggregate_master[shared_rows_aggregate_master[self.ALT_DEPTH] == 0])

        # the number of rows we should impute across all samples and generate a refDepth for
        n_imputed_zero_altDepth = n_zero_altDepth_aggregate - (n_zero_altDepth_from_calls + n_zero_altDepth_from_master)

        # number of total rows we should have in aggregate
        n_unique_chr_pos = len(self.unique_chr_pos) * len(self.populations)


        # the number of values we should impute are the total number of <chromosome><position> pairs
        # minus the number of row we found in the master df, and calls df, while accounting for duplicates
        n_missing_entry_master_or_calls = n_unique_chr_pos - len(self.master_df)

        # rows pulled from calls
        rows_pulled_from_calls = pd.concat([shared_rows_aggregate_calls, self.master_df]).drop_duplicates(subset=on_list_no_alt, keep=False)


        if self.metrics_file:

            pdf = MPN_AML_METRICS_PDF(filename=self.metrics_file,
                                      title="MPN-AML-Aggregator Metrics",
                                      details="Master xls: %s\n Calls xls: %s\n Populations xls: %s\n Aggregated xls: %s" \
                                      % (master_xls, calls_xls, populations_xls, aggregated_xls))


            pdf.add_details(
                conditions=[

                    all(self.aggregated_df.groupby(self.CHR_POS)[self.GENE].nunique().eq(1)) and \
                    (self.aggregated_df.groupby(self.CHR_POS)[self.GENE].value_counts().index[0] == self.master_df.groupby(self.CHR_POS)[self.GENE].value_counts().index[0]),

                    len(self.aggregated_df) == len(self.master_df[self.CHR_POS].unique()) * len(self.populations),

                    len(shared_rows_aggregate_master) == len(self.master_df),

                    n_imputed_zero_altDepth == n_missing_entry_master_or_calls,

                    True,

                    len(imputed_rows) == (n_unique_chr_pos-len(shared_rows_aggregate_calls)),

                    len(rows_pulled_from_calls) ==  (n_unique_chr_pos - len(imputed_rows) - len(self.master_df))

                ],
                statements=[

                    "all unique <chromosome><position> pairs have the same gene in both the aggregated xls and master xls",

                    "number of rows in aggregate / ((# of samples) * (# of unique <chromosome><position> pairs)) = %d/%d"
                         % (len(self.aggregated_df), len(self.master_df[self.CHR_POS].unique()) * len(self.populations)),

                    "number of matching rows between aggregate and master / number of rows in master = %d/%d"
                         % (len(shared_rows_aggregate_master), len(self.master_df)),

                    "number of <chromosome><position> pairs with altDepth = 0 / remaining unmatched rows (from either calls or master) = %d/%d"
                         % (n_imputed_zero_altDepth, n_missing_entry_master_or_calls),

                    "number of matching rows between aggregate and calls / number of rows in calls = %d/%d"
                          % (len(shared_rows_aggregate_calls), len(self.calls_df)),

                    "number of rows not found in either calls or master xls (have imputed refDepth) / total number of rows = %d/%d" % (len(imputed_rows), n_unique_chr_pos),

                    "number of rows pulled from calls df / total number of rows = %d/%d" % (len(rows_pulled_from_calls), n_unique_chr_pos)

                ],
                title="Aggregation Details"

            )


            if len(imputed_rows) != 0:

                pdf.add_plot(
                    x_data=imputed_rows.drop_duplicates(self.CHR_POS)[self.CHR_POS], y_data=imputed_rows.drop_duplicates(self.CHR_POS)[self.REF_DEPTH],
                    suptitle="Imputed Ref Depth Per Chromosome-Position", title="Rows Missing From Master xls",
                    xlabel="Chromosome_Position", ylabel="Imputed Ref Depth",
                    xtick_rot=90
                )

                # pdf.add_table(
                #     data=imputed_rows[table_columns],
                #     title="Rows with imputed refDepth",
                #     tight_layout=False
                # )


            # plot per sample metrics
            pbar = tqdm(range(0, len(self.populations)))

            for idx in pbar:

                pop = self.populations[idx]

                pbar.set_description("Generating metrics for %s" % pop)

                # dataframes used to generate plots
                sample_agg_df = self.aggregated_df[self.aggregated_columns].loc[self.aggregated_df[self.SAMPLE_NAMES].eq(pop)]
                sample_shared_agg_calls_w_alt = shared_rows_aggregate_calls.loc[shared_rows_aggregate_calls[self.SAMPLE_NAMES].eq(pop)]
                sample_master_df = self.master_df.loc[self.master_df[self.SAMPLE_NAMES].eq(pop)]
                sample_imputed_rows = pd.concat([sample_agg_df, sample_shared_agg_calls_w_alt]).drop_duplicates(subset=on_list_no_alt, keep=False)

                # table of rows that are in calls xls, but not in master xls
                rows_of_sample_not_in_master = pd.concat([sample_master_df, sample_shared_agg_calls_w_alt]).drop_duplicates(subset=[self.CHR_POS], keep=False)


                pdf.add_plot(
                    x_data=sample_agg_df[self.CHR_POS], y_data=sample_agg_df[self.VAF],
                    suptitle="VAF per Chromosome-Position (Aggregate of master,calls,imputed)", title="Sample %s" % pop,
                    xlabel="Chromosome_Position", ylabel="Variant Allele Frequency (VAF)",
                    caption="# of aggregated Chromosome_Position: %d\n # from master xls: %d\n # from calls xls: %d\n # imputed: %d" \
                                % (len(sample_agg_df), len(sample_master_df), len(sample_shared_agg_calls_w_alt)-len(sample_master_df), len(sample_imputed_rows)),
                    xtick_rot=90, ylim=[0.0, 1.0]
                )

                if len(rows_of_sample_not_in_master) != 0:

                    pdf.add_plot(
                        x_data=rows_of_sample_not_in_master[self.CHR_POS], y_data=rows_of_sample_not_in_master[self.VAF],
                        suptitle="VAF per Chromosome-Position (from calls xls)", title="Sample %s" % pop,
                        xlabel="Chromosome_Position", ylabel="Variant Allele Frequency (VAF)",
                        caption="# from calls xls: %d" % len(rows_of_sample_not_in_master),
                        xtick_rot=90, ylim=[0.0, 1.0]
                    )

                    pdf.add_table(
                        data=rows_of_sample_not_in_master[table_columns],
                        title="Entries from calls xls for sample %s" % pop
                    )

            pdf.render()
