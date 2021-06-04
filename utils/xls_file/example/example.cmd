# Example of a test command file for generate an aggregated xlsx file
# Usage (assuming you will read and write files in the data directory):
# python run_aggregation.py -f example/example.cmd -i $DATA_DIR/ -j $DATA_DIR/



-m MATS08_sp.xlsx SNVs
-c MATS08_sq.xlsx Sheet1
-p MATS08_sp.xlsx Populations None
-o agg_testout.xlsx Sheet1
-a MPN_AML_Aggregator
