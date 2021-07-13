import operator
import argparse

from modify_pop import load_excel, load_csv, save_excel, select_pops, remove_pops

# to run an example, use the following command:
#   python3 $UTILS_DIR/ssm_file/utils/run_modify_pop.py -i example.output.ssm -o example.modified.ssm -d $DATA_DIR/example/results/ -a \> 0.5 -m RM_VARS_BY_VAF


MOD_METHODS = {
    "SELECT_POPS": select_pops,
    "REMOVE_POPS": remove_pops
}


def _parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(

        description='Update an ssm file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter

    )

    parser.add_argument('-p', '--population-file', nargs='+', help='Population file <file_name> <sheet_name> <header>', required=True)
    parser.add_argument('-o', '--out-file', nargs='+', help='Modified population file <file_name> <sheet_name> <header>', required=True)
    parser.add_argument('-c', '--mod-file', help='Csv of populations to modify.')
    parser.add_argument('-m', '--mod-method', help='Modification method to be applied to ssm file.', choices=tuple(MOD_METHODS.keys()), required=True)


    args = parser.parse_args()

    return args


def main():
    """
    Performs checks on command line arguments, then attempts to process file.
    """
    args = _parse_args()

    if args.population_file[-1] == 'None':
        args.population_file[-1] = None

    if args.out_file[-1] == 'None':
        args.out_file[-1] = None

    # apply method
    save_excel(

        MOD_METHODS[args.mod_method](
            load_excel(*args.population_file),
            load_csv(args.mod_file)
        ),
        *args.out_file

    )


if __name__ == '__main__':
  main()
