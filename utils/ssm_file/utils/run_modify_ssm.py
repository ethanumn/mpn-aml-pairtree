import operator
import argparse

from modify_ssm import load_ssm, load_csv, save_ssm, remove_vars_by_vaf, organize_vars_by_vaf, scale_counts, separate_garbage, keep_vars_by_name

# to run an example, use the following command:
#   python3 $UTILS_DIR/ssm_file/utils/run_modify_ssm.py -i example.output.ssm -o example.modified.ssm -d $DATA_DIR/example/results/ -a \> 0.5 -m RM_VARS_BY_VAF


MOD_METHODS = {
    "RM_VARS_BY_VAF": remove_vars_by_vaf,
    "ORG_VARS_BY_VAF": organize_vars_by_vaf,
    "SCALE_COUNTS": scale_counts,
    "SEPARATE_GARBAGE": separate_garbage,
    "KEEP_VARS_BY_NAME": keep_vars_by_name
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
    parser.add_argument('-p', '--params-file', help='Name of params file to reference.', default=None)
    parser.add_argument('-d', '--directory', help='Directory to read/write files from.')
    parser.add_argument('-a', '--args', nargs='+', help='Additional arguments to pass to modification method.')
    parser.add_argument('-m', '--mod-method', help='Modification method to be applied to ssm file.', choices=tuple(MOD_METHODS.keys()), required=True)
    parser.add_argument('-n', '--names-fn', help='File containing names to keep')


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

    # extract args passed

    # keep names
    if args.names_fn and args.mod_method == "KEEP_VARS_BY_NAME":
        args.args = [open(args.names_fn).read().splitlines()]
    # for separate garbage
    elif len(args.args) == 0 and args.params_file:
        args.args = args.params_file

    # for scale counts
    elif len(args.args) == 1:
        if args.params_file:
            args.args = [args.params_file, args.args[0]]
        else:
            args.args = [int(args.args[0])]

    # for rm vaf
    elif len(args.args) == 2:
        args.args = [OPERATORS[args.args[0]], float(args.args[1])]

    # for organize vaf
    elif len(args.args) == 3:
        args.args = [OPERATORS[args.args[0]], float(args.args[1]), args.args[2]]

    # apply method
    save_ssm(

        MOD_METHODS[args.mod_method](load_ssm(args.in_file), *args.args),
        args.out_file

    )


if __name__ == '__main__':
  main()
