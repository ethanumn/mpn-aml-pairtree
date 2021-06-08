# Example of a test command file for generate .ssm files
# Usage (assuming your in-files and out-files will go in the data directory):
# python3 run_processor.py -f example/example.cmd -d $DATA_DIR/example/


-i example.aggregated.xlsx
-o example.output.ssm
-p MPN_AML_Processor
