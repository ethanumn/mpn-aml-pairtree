import pandas as pd

class SSM_Base_Processor:
    """
    Base class that can be used to generate ssm files.

    For an example of its use, see 'create_aml_ssm_file.py'.
    """

    def __init__(self, in_file="", out_file="", write_out_file=True, write_out_params=True):

        # set up everything necessary to read/process/write
        self._init_constants()
        self._init_variables()

        # start processing if we have all of the information we need (I/O file names)
        if in_file:

            # read and set the in-file
            self.read_in_file(in_file)

            # run all processing functions
            self.process()

            # write SSM file
            if write_out_file and out_file:
                self.write_out_file(out_file)

            # write the params file
            if write_out_params and out_file:
                self.write_out_params(out_file.replace(".ssm", ".params.json"))


    def _init_constants(self):

        # .xlsx aggregated file columns
        self.CHR = "Chr"
        self.POSITION = "Position"
        self.CHR_POS = "chr_pos"
        self.REF_DEPTH = "refDepth"
        self.ALT_DEPTH = "altDepth"
        self.SAMPLE_NAMES = "sampleNames"
        self.VAF = "VAF"

        # .json columns
        self.SAMPLES = "samples"
        self.CLUSTERS = "clusters"
        self.GARBAGE = "garbage"

        # .ssm column
        self.P_SIGNATURE = "p_.*"

        self.COL_NAME = "name"

        self.COL_ID = "id"

        self.COL_VAR_READS = "var_reads"

        self.COL_TOTAL_READS = "total_reads"

        self.COL_VAR_READ_PROB = "var_read_prob"

        self.COL_ORDER = [

            self.COL_ID,
            self.COL_NAME,
            self.COL_VAR_READS,
            self.COL_TOTAL_READS,
            self.COL_VAR_READ_PROB

        ]


    def _init_variables(self):

        # initializations
        self.in_file = None
        self.out_file = None

        self.in_df = None
        self.out_df = pd.DataFrame()

        self.processing_functions = [
            # all functions used to translate input file to SSM file
        ]


    def read_in_file(self, in_file=""):
        """
        Can be used to either read the in-file, or both set and read the in-file.
        This method may want to be used outside of t
        """

        if in_file:

            self.in_file = in_file

        if self.in_file:

            self.in_df = pd.read_excel(self.in_file)


    def write_out_file(self, out_file=""):

        if out_file:

            self.out_file = out_file

        if self.out_file and isinstance(self.out_df, pd.DataFrame):

            if set(self.COL_ORDER).issubset(self.out_df.columns):

                self.out_df[self.COL_ORDER].to_csv(self.out_file, sep="\t", index=False)


    def write_out_params(self, params_file=""):

        import json

        def convert(obj):
            """
            Overrides default function when something isn't json serializable
            """
            import numpy as np

            if isinstance(obj, np.ndarray): return obj.tolist()
            else: raise TypeError

        if params_file:

            self.params_file = params_file

        if self.params_file and isinstance(self.out_df, pd.DataFrame):

            if set(self.COL_ORDER).issubset(self.out_df.columns):

                # regroups dataframe by sample name, then collects all unique ids per sample
                params_dict = self.in_df.groupby(self.SAMPLE_NAMES).apply(lambda sample_name: sample_name[self.COL_ID].unique())

                with open(params_file, 'w') as outfile:

                    json_dict = json.dumps(
                                    {
                                      self.SAMPLES: [params_dict.keys().to_numpy()],
                                      self.CLUSTERS: [params_dict.values],
                                      self.GARBAGE: self.garbage_mutations()
                                    },
                                    default=convert
                                )

                    outfile.write(json_dict)


    def garbage_mutations(self):
        """
        This may be replaced or implemented at the subclass level to detect garbage mutations to ignore.
        """
        return []


    def process(self):
        """
        Collects all processing functions (functions that match r'p\_.*'), calls all processing
        functions, and then writes the resulting pandas dataframe to a file
        """

        from tqdm import tqdm

        # collect all processing functions
        self._aggregate_processing_functions()

        # set up progress bar with all processing functions
        pbar = tqdm(self.processing_functions)

        # call all processing functions
        for function in pbar:

            pbar.set_description("Running %s" % function.__name__.lstrip(self.P_SIGNATURE))

            function()

            pbar.set_description("Completed.")

        pbar.set_description("All processing complete.")


    def _aggregate_processing_functions(self):
        """
        Searches through object attributes and methods for functions that match r'p\_.*'
        """

        import re

        regex = re.compile(self.P_SIGNATURE)

        # finds all attributes which match our prefix signature (should be processing functions)
        processing_functions = list(filter(regex.match, dir(self)))

        self.processing_functions = [getattr(self, f_name) for f_name in processing_functions]
