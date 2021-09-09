import operator
import argparse

from modify_tsv import load_ssm, load_csv, save_ssm, match_tsv_to_ssm

# to run an example, use the following command:
#  python3 $UTILS_DIR/ssm_file/utils/tsv_file/tsv_to_ssm.py -s example.ssm -t example.tsv -p example.params.json -m MATCH_TSV_TO_SSM


MOD_METHODS = {
    "MATCH_TSV_TO_SSM": match_tsv_to_ssm,
}


def _parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(

        description='Update an ssm file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-s', '--ssm-file', help='Name of ssm file to modify.', required=True)
    parser.add_argument('-t', '--tsv-file', help='Name of tsv file to modify.', required=True)
    parser.add_argument('-p', '--params-file', help='Name of params file to modify.', default=None)
    parser.add_argument('-d', '--directory', help='Directory to read/write files from.')
    parser.add_argument('-m', '--mod-method', help='Modification method to be applied to ssm file.', choices=tuple(MOD_METHODS.keys()), required=True)


    args = parser.parse_args()

    return args


def main():
    """
    Performs checks on command line arguments, then attempts to process all files.
    """
    args = _parse_args()


    # if we want to pass in params file argument
    if args.params_file == "None":
        args.params_file = None

    # append directory to in_file/out_file if the argument was passed
    if args.directory:
        args.in_file = args.directory + args.in_file
        args.out_file = args.directory + args.out_file

    if not args.params_file or not args.ssm_file or not args.tsv_file:
        raise FileNotFoundError("Did not pass in params, ssm or tsv file")

    MOD_METHODS[args.mod_method](args.tsv_file, args.ssm_file, args.params_file),




if __name__ == '__main__':
  main()
