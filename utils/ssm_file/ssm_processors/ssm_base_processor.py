import pandas as pd
import sys, os

sys.path.append(os.environ["UTILS_DIR"] + "/common")

from mpn_aml_columns import *
from ssm_columns import *

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

        # .ssm column
        self.P_SIGNATURE = "p_.*"

        # column order to output
        self.COL_ORDER = [

            COL_ID,
            COL_NAME,
            COL_VAR_READS,
            COL_TOTAL_READS,
            COL_VAR_READ_PROB

        ]


    def _init_variables(self):

        # initializations
        self.in_file = None
        self.out_file = None

        self.in_df = None
        self.processed_df = pd.DataFrame()
        self.out_df = None

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


    def format_out_df(self):
        """
        Empty base function that should be overriden in child class that will translate the processed_df
        into a dataframe which will be written out as a .ssm file
        """
        pass


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

        if self.params_file and isinstance(self.in_df, pd.DataFrame):

            if set([SAMPLE_NAMES]).issubset(self.in_df.columns):

                with open(self.params_file, 'w') as outfile:

                    json_dict = json.dumps(
                                    {
                                      SAMPLES: [self.in_df[SAMPLE_NAMES].unique()],
#                                      CLUSTERS: [],
                                      GARBAGE: self.garbage_mutations()
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
        functions, then transforms the processed dataframe into a dataframe that can written as a .ssm
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

        # translate processed_df into out_df for writing to file
        self.format_out_df()


    def _aggregate_processing_functions(self):
        """
        Searches through object attributes and methods for functions that match r'p\_.*'
        """

        import re

        regex = re.compile(self.P_SIGNATURE)

        # finds all attributes which match our prefix signature (should be processing functions)
        processing_functions = list(filter(regex.match, dir(self)))

        self.processing_functions = [getattr(self, f_name) for f_name in processing_functions]
