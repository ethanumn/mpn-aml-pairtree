import operator
import argparse

from modify_ssm import load_ssm, save_ssm, remove_vars_by_vaf, organize_vars_by_vaf, scale_counts

# to run an example, use the following command:
#   python3 $UTILS_DIR/ssm_file/utils/run_modify_ssm.py -i example.output.ssm -o example.modified.ssm -d $DATA_DIR/example/results/ -a \> 0.5 -m RM_VARS_BY_VAF


MOD_METHODS = {
    "RM_VARS_BY_VAF": remove_vars_by_vaf,
    "ORG_VARS_BY_VAF": organize_vars_by_vaf,
    "SCALE_COUNTS": scale_counts
}

OPERATORS = {
    "==": operator.eq,
    "<" : operator.lt,
    ">" : operator.gt,
    ">=": operator.ge,
    "<=": operator.le,
    "!=": operator.ne
}


def _parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(

        description='Update an ssm file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-i', '--in-file', help='Name of ssm file to modify.', required=True)
    parser.add_argument('-o', '--out-file', help='Name of modified ssm file.', required=True)
    parser.add_argument('-d', '--directory', help='Directory to read/write files from.')
    parser.add_argument('-a', '--args', nargs='+', help='Additional arguments to pass to modification method.')
    parser.add_argument('-m', '--mod-method', help='Modification method to be applied to ssm file.', choices=tuple(MOD_METHODS.keys()), required=True)


    args = parser.parse_args()

    return args


def main():
    """
    Performs checks on command line arguments, then attempts to process all files.
    """
    args = _parse_args()

    # append directory to in_file/out_file if the argument was passed
    if args.directory:
        args.in_file = args.directory + args.in_file
        args.out_file = args.directory + args.out_file


    # extract args passed
    if len(args.args) == 1:
        args.args = [int(args.args[0])]

    elif len(args.args) == 2:
        args.args = [OPERATORS[args.args[0]], float(args.args[1])]

    elif len(args.args) == 3:
        args.args = [OPERATORS[args.args[0]], float(args.args[1]), args.args[2]]

    # apply method
    save_ssm(

        MOD_METHODS[args.mod_method](load_ssm(args.in_file), *args.args),
        args.out_file

    )


if __name__ == '__main__':
  main()
