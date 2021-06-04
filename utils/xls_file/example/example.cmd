cd# Example of a test command file for generate an aggregated xlsx file
# Usage (assuming you will read and write files in the data directory):
# python run_aggregator.py -f example/example.cmd -i $DATA_DIR/example/ -j $DATA_DIR/example/


-m example.master.xlsx Sheet1
-c example.calls.xlsx Sheet1
-p example.populations.xlsx Sheet1 None
-o example.aggregated.xlsx Sheet1
-a MPN_AML_Aggregator
