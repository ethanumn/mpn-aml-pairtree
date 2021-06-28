import argparse
import os
import sys

sys.path.append(os.environ["UTILS_DIR"] + "/ssm_file/ssm_processors")

from ssm_base_processor import SSM_Base_Processor
from mpn_aml_processor import MPN_AML_Processor


# NEED to add any processor you might want to use
PROCESSORS = {
    SSM_Base_Processor.__name__ : SSM_Base_Processor,
    MPN_AML_Processor.__name__ : MPN_AML_Processor
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


def _parse_args(processor_choices=[]):
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(

        description='Process a custom file type into an .ssm file readable by Pairtree',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-f', '-arg-file', type=open, action=LoadFromFile)

    parser.add_argument('-i', '--in-files', nargs='+', help='List of in-files to process')
    parser.add_argument('-o', '--out-files', nargs='+', help='List of out-file names (must occur in same order as corresponding in-file)')
    parser.add_argument('-p', '--processors', nargs='+', help='List of processor types for each corresponding in-file', choices=tuple(processor_choices))
    parser.add_argument('-d', '--directories', nargs='+', help='List of directories to read/write files from')

    args = parser.parse_args()

    return args


def run_processors(processors, in_files, out_files, directories):
    """
    Runs all processors dependent on what arguments are passed via the command line
    """
    if (processors == None) or (in_files == None) or (out_files) == None:
        raise argparse.ArgumentTypeError('did not pass in one or more arguments (--in-files, --out-files, --processors)')


    # number of in-files needs to match the number of out-files
    if len(in_files) != len(out_files):
        raise argparse.ArgumentTypeError('in-file count does not match out-file count')

    # concatenate directories with file names if necessary
    if directories != None:

        if len(directories) == 1:
            for idx in range(0, len(in_files)):
                in_files[idx] = directories[0] + in_files[idx]
                out_files[idx] = directories[0] + out_files[idx]

        elif len(directories) == len(in_files):
            for idx in range(0, len(in_files)):
                in_files[idx] = directories[idx] + in_files[idx]
                out_files[idx] = directories[idx] + out_files[idx]
        else:
            raise argparse.ArgumentTypeError('in-file count does not match directories count')



    # if we only have one processor, use it for all of our files
    if len(processors) == 1:
        for in_file, out_file in zip(in_files, out_files):
            processors[0](in_file, out_file)

    else:
        if len(in_files) != len(processors):
            raise argparse.ArgumentTypeError('in-file count does not match processor count')

        for processor, in_file, out_file in zip(processors, in_files, out_files):
            processor(in_file, out_file)


def main():
    """
    Performs checks on command line arguments, then attempts to process all files.
    """
    args = _parse_args(PROCESSORS.keys())

    run_processors([PROCESSORS[processor_name] for processor_name in args.processors],
                   args.in_files,
                   args.out_files,
                   args.directories)


if __name__ == '__main__':
  main()
