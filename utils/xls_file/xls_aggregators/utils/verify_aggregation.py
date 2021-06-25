import pandas as pd
import sys, os

sys.path.append(os.environ["UTILS_DIR"] + "/pdf_templates")
sys.path.append(os.environ["UTILS_DIR"] + "/common")


from mpn_aml_metrics_pdf import MPN_AML_METRICS_PDF
from mpn_aml_columns import *

def verify_aggregation(metrics_file,
                       aggregated_df,
                       primary_df,
                       calls_df,
                       populations,
                       unique_chr_pos,
                       aggregated_columns,
                       primary_xls="",
                       calls_xls="",
                       populations_xls="",
                       aggregated_xls=""):
    """
    Verify aggregation by producing a pdf that contains the following:
        - Verification of the number of rows in the different dataframes and some unique corner cases
        - Plot/Table of imputed refDepth values for any unique <chromosome><position>
          that should have been in the primary xls but was not
        - Plot of VAF of each <chromosome><position> pair per sample
        - Plot/Table of each <chromosome><position> pair pulled from the calls xls per sample
    """
    from tqdm import tqdm


    # list of columns to join on - do not use VAF or any floating point value since the changes in precision
    # do not allow for a proper join. DO NOT use altDepth when comparing between calls_df and aggregate_df because we've zero'd out
    # the altDepth of those rows we pulled from calls_df into the aggregate_df
    on_list = [CHR_POS, CHR, POSITION, SAMPLE_NAMES, ALT_DEPTH, REF_DEPTH]
    on_list_no_alt = [CHR_POS, CHR, POSITION, SAMPLE_NAMES, REF_DEPTH]
    table_columns = [SAMPLE_NAMES, CHR_POS, POSITION, ALT_DEPTH, REF_DEPTH]

    # inner join between primary_df and aggregate_df
    shared_rows_aggregate_primary = pd.merge(primary_df, aggregated_df, how="inner", on=on_list, suffixes=("_x", ""))

    # inner join between calls_df and primary_df
    shared_rows_primary_calls = pd.merge(calls_df.drop(columns=[VAF]), primary_df, how="inner", on=on_list, suffixes=("_x", ""))

    # inner join between aggregate_df and calls_df
    shared_rows_aggregate_calls = pd.merge(aggregated_df, calls_df.drop(columns=[VAF]), how="inner",
                                                on=on_list_no_alt, suffixes=("_x", "")).dropna()

    # rows that have ref depth imputed
    imputed_rows = pd.concat([aggregated_df, shared_rows_aggregate_calls]).drop_duplicates(subset=on_list_no_alt, keep=False)

    # number of possible imputed rows
    n_zero_altDepth_aggregate = len(aggregated_df[aggregated_df[ALT_DEPTH] == 0])

    # number of rows that came from calls dataframe that have an altDepth = 0
    n_zero_altDepth_from_calls = len(shared_rows_aggregate_calls[shared_rows_aggregate_calls[ALT_DEPTH] == 0])

    # number of rows that came from primary dataframe that have an altDepth = 0
    n_zero_altDepth_from_primary = len(shared_rows_aggregate_primary[shared_rows_aggregate_primary[ALT_DEPTH] == 0])

    # the number of rows we should impute across all samples and generate a refDepth for
    n_imputed_zero_altDepth = n_zero_altDepth_aggregate - (n_zero_altDepth_from_calls + n_zero_altDepth_from_primary)

    # number of total rows we should have in aggregate
    n_unique_chr_pos = len(unique_chr_pos) * len(populations)


    # the number of values we should impute are the total number of <chromosome><position> pairs
    # minus the number of row we found in the primary df, and calls df, while accounting for duplicates
    n_missing_entry_primary_or_calls = n_unique_chr_pos - len(primary_df)

    # rows pulled from calls
    rows_pulled_from_calls = pd.concat([shared_rows_aggregate_calls, primary_df]).drop_duplicates(subset=on_list_no_alt, keep=False)

    rows_pulled_from_calls[VAF] = rows_pulled_from_calls[ALT_DEPTH] / (rows_pulled_from_calls[ALT_DEPTH] + rows_pulled_from_calls[REF_DEPTH])

    if metrics_file:

        pdf = MPN_AML_METRICS_PDF(filename=metrics_file,
                                  title="MPN-AML-Aggregator Metrics",
                                  details="Primary xls: %s\n Calls xls: %s\n Populations xls: %s\n Aggregated xls: %s" \
                                  % (primary_xls, calls_xls, populations_xls, aggregated_xls))


        # add some debug statements to verify our aggregation worked properly
        pdf.add_details(
            conditions=[

                all(aggregated_df.groupby(CHR_POS)[GENE].nunique().eq(1)) and \
                (aggregated_df.groupby(CHR_POS)[GENE].value_counts().index[0] == primary_df.groupby(CHR_POS)[GENE].value_counts().index[0]),

                len(aggregated_df) == len(primary_df[CHR_POS].unique()) * len(populations),

                len(shared_rows_aggregate_primary) == len(primary_df),

                n_imputed_zero_altDepth == n_missing_entry_primary_or_calls,

                True,

                len(imputed_rows) == (n_unique_chr_pos-len(shared_rows_aggregate_calls)),

                len(rows_pulled_from_calls) ==  (n_unique_chr_pos - len(imputed_rows) - len(primary_df)),

                True

            ],
            statements=[

                "all unique <chromosome><position> pairs have the same gene in both the aggregated xls and primary xls",

                "number of rows in aggregate / ((# of samples) * (# of unique <chromosome><position> pairs)) = %d/%d"
                     % (len(aggregated_df), len(primary_df[CHR_POS].unique()) * len(populations)),

                "number of matching rows between aggregate and primary / number of rows in primary = %d/%d"
                     % (len(shared_rows_aggregate_primary), len(primary_df)),

                "number of <chromosome><position> pairs with altDepth = 0 / remaining unmatched rows (from either calls or primary) = %d/%d"
                     % (n_imputed_zero_altDepth, n_missing_entry_primary_or_calls),

                "number of matching rows between aggregate and calls / number of rows in calls = %d/%d"
                      % (len(shared_rows_aggregate_calls), len(calls_df)),

                "number of rows not found in either calls or primary xls (have imputed refDepth) / total number of rows = %d/%d" % (len(imputed_rows), n_unique_chr_pos),

                "number of rows pulled from calls df / total number of rows = %d/%d" % (len(rows_pulled_from_calls), n_unique_chr_pos),

                "number of unique chromosome_position where a refDepth imputation had to be done / number of unique chromosome_position = %d/%d" % (len(imputed_rows[CHR_POS].unique()), len(primary_df[CHR_POS].unique()))
            ],
            title="Aggregation Details"

        )


        # plot refDepth of imputed rows
        if len(imputed_rows) != 0:

            pdf.add_plot(
                x_data=imputed_rows.drop_duplicates(CHR_POS)[CHR_POS], y_data=imputed_rows.drop_duplicates(CHR_POS)[REF_DEPTH],
                suptitle="Imputed Ref Depth Per Chromosome-Position", title="Rows Missing From primary xls",
                xlabel="Chromosome_Position", ylabel="Imputed Ref Depth",
                xtick_rot=90
            )

            # pdf.add_table(
            #     data=imputed_rows[table_columns],
            #     title="Rows with imputed refDepth",
            #     tight_layout=False
            # )

        # table of variants with VAF > 0.5 from primary xls
        pdf.add_table(
            data=primary_df.loc[primary_df[VAF] > 0.5, [SAMPLE_NAMES, CHR, POSITION, GENE, VAF, ALT_DEPTH, REF_DEPTH]],
            title="Variants from filtered xls with a VAFs > 0.5"
        )


        vaf_threshold = .1
        pdf.add_table(
            data=rows_pulled_from_calls.loc[rows_pulled_from_calls[VAF] > vaf_threshold, [SAMPLE_NAMES, CHR, POSITION, GENE, VAF, ALT_DEPTH, REF_DEPTH]],
            title=("Variants from calls xls with a VAF > %.2f" % vaf_threshold)
        )


        # initialize sample overview dataframe {"sample": ?, "variants_in_sample": ?, "variants_from_primary": ?, "variants_from_calls": ?, variants_imputed: ?}
        sample_overview_df = pd.DataFrame()

        # plot per sample metrics
        pbar = tqdm(range(0, len(populations)))

        for idx in pbar:

            pop = populations[idx]

            pbar.set_description("Generating metrics for %s" % pop)

            # dataframes used to generate plots
            sample_agg_df = aggregated_df[aggregated_columns].loc[aggregated_df[SAMPLE_NAMES].eq(pop)]
            sample_shared_agg_calls_w_alt = shared_rows_aggregate_calls.loc[shared_rows_aggregate_calls[SAMPLE_NAMES].eq(pop)]
            sample_primary_df = primary_df.loc[primary_df[SAMPLE_NAMES].eq(pop)]
            sample_imputed_rows = pd.concat([sample_agg_df, sample_shared_agg_calls_w_alt]).drop_duplicates(subset=on_list_no_alt, keep=False)

            # table of rows that are in calls xls, but not in primary xls
            rows_of_sample_not_in_primary = pd.concat([sample_primary_df, sample_shared_agg_calls_w_alt]).drop_duplicates(subset=CHR_POS, keep=False)

            # append row to sample_overview_df which gives the count of where information was pulled from for each sample
            sample_overview_rows = [SAMPLE_NAMES, "Total Variants in Sample", "Variants From Filtered", "Variants From Calls", "Variants Imputed" ]
            sample_overview_df = sample_overview_df.append(
                dict(zip(
                        sample_overview_rows,
                        [pop, len(sample_agg_df), len(sample_primary_df), len(rows_of_sample_not_in_primary), len(sample_imputed_rows)]
                )), ignore_index=True
            )

            # plot of all VAFs per sample
            pdf.add_plot(
                x_data=sample_agg_df[CHR_POS], y_data=sample_agg_df[VAF],
                suptitle="VAF per Chromosome-Position (Aggregate of primary,calls,imputed)", title="Sample %s" % pop,
                xlabel="Chromosome_Position", ylabel="Variant Allele Frequency (VAF)",
                caption="# of aggregated Chromosome_Position: %d\n # from primary xls: %d\n # from calls xls: %d\n # imputed: %d" \
                            % (len(sample_agg_df), len(sample_primary_df), len(sample_shared_agg_calls_w_alt)-len(sample_primary_df), len(sample_imputed_rows)),
                xtick_rot=90, ylim=[0.0, 1.0]
            )


            if len(rows_of_sample_not_in_primary) != 0:

                # plot of VAFs from calls xls per sample
                pdf.add_plot(
                    x_data=rows_of_sample_not_in_primary[CHR_POS],
                    y_data=rows_of_sample_not_in_primary[ALT_DEPTH] / (rows_of_sample_not_in_primary[ALT_DEPTH]+rows_of_sample_not_in_primary[REF_DEPTH]),
                    suptitle="VAF per Chromosome-Position (from calls xls)", title="Sample %s" % pop,
                    xlabel="Chromosome_Position", ylabel="Variant Allele Frequency (VAF)",
                    caption="# from calls xls: %d" % len(rows_of_sample_not_in_primary),
                    xtick_rot=90, ylim=[0.0, 1.0]
                )

                # table of VAFs from calls xls per sample
                pdf.add_table(
                    data=rows_of_sample_not_in_primary[table_columns],
                    title="Entries from calls xls for sample %s" % pop
                )

            if len(sample_imputed_rows) != 0:
                # table of imputed row per sample
                pdf.add_table(
                    data=sample_imputed_rows[table_columns],
                    title="Imputed rows for sample %s" % pop
                )



        if len(sample_overview_df) != 0:

            # plot table of overview of which xls variants came from for each sample
            sample_overview_df.loc[:, sample_overview_df.columns != SAMPLE_NAMES] = sample_overview_df.loc[:, sample_overview_df.columns != SAMPLE_NAMES].astype("int64")
            pdf.add_table(data=sample_overview_df[sample_overview_rows],
                          title="Source of variant data for each sample")


        pdf.render()
