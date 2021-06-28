import argparse
import os
import sys

sys.path.append(os.environ["UTILS_DIR"] + "/xls_file/xls_aggregators")

from mpn_aml_aggregator import MPN_AML_Aggregator, IMPUTE_AVG, IMPUTE_ZERO


# NEED to add any aggregator you might want to use
aggregator_dict = {
    MPN_AML_Aggregator.__name__ : MPN_AML_Aggregator
}


class LoadFromFile(argparse.Action):
    """
    Loads command line arguments from file.
    """
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:

            args = []

            # we want to ignore all lines that start with a #
            for line in f.readlines():
                if not line.startswith("#"):
                    args += line.split()

            # parse arguments in the file and store them in the target namespace
            parser.parse_args(args, namespace)


def _parse_args(aggregator_choices=[]):
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(

        description='Process a custom file type into an .ssm file readable by Pairtree',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-f', '-arg-file', type=open, action=LoadFromFile)

    parser.add_argument('-m', '--primary-file', nargs='+', help='primary file <file_name> <sheet_name>')
    parser.add_argument('-c', '--call-file', nargs='+', help='Call file <file_name> <sheet_name>')
    parser.add_argument('-o', '--output-file', nargs='+', help='Output file <file_name> <sheet_name>')
    parser.add_argument('-p', '--population-file', nargs='+', help='Population file <file_name> <sheet_name> <header>')
    parser.add_argument('-d', '--metrics-file', default="", help='File to output aggregation metrics to <file_name>')
    parser.add_argument('-a', '--aggregator', help='Aggregator to run on files', choices=tuple(aggregator_choices))
    parser.add_argument('-i', '--input-directory', help='Directory to read primary/call files from')
    parser.add_argument('-j', '--output-directory', help='Directory to write aggregated file to')
    parser.add_argument('-t', '--impute-technique', default=IMPUTE_ZERO, help='Technique to use for imputing missing values', choices=(IMPUTE_AVG, IMPUTE_ZERO))
    args = parser.parse_args()

    return args


def run_aggregators(aggregator, primary_file, call_file, population_file, output_file, metrics_file, input_directory, output_directory, impute_technique):
    """
    Runs all aggregators dependent on what arguments are passed via the command line
    """

    # check to make sure we've been passed files
    if len(primary_file) == 0:
        raise argparse.ArgumentTypeError('No primary file specified')

    if len(call_file) == 0:
        raise argparse.ArgumentTypeError('No call file specified')

    if len(output_file) == 0:
        raise argparse.ArgumentTypeError('No output file specified')

    if len(population_file) == 0:
        raise argparse.ArgumentTypeError('No population file specified')


    # concatenate input directories with file names if necessary
    if input_directory != None:
        primary_file[0] = input_directory + primary_file[0]
        call_file[0] = input_directory + call_file[0]
        population_file[0] = input_directory + population_file[0]


    # concatenate output directory with file names if necessary
    if output_directory != None:
        output_file[0] = output_directory + output_file[0]
        metrics_file = output_directory + metrics_file



    # workaround for passing header
    if population_file[-1] == 'None':
        population_file[-1] = None

    if primary_file[-1] == 'None':
        primary_file[-1] = None

    if call_file[-1] == 'None':
        call_file[-1] = None

    # if we only have one aggregator, use it for all of our files
    if aggregator != None:
        aggregator(primary_file, call_file, population_file, output_file, metrics_file, impute_technique)



def main():
    """
    Performs checks on command line arguments, then attempts to process all files.
    """
    args = _parse_args(aggregator_dict.keys())

    run_aggregators(aggregator_dict[args.aggregator],
                    args.primary_file,
                    args.call_file,
                    args.population_file,
                    args.output_file,
                    args.metrics_file,
                    args.input_directory,
                    args.output_directory,
                    args.impute_technique)


if __name__ == '__main__':
  main()
