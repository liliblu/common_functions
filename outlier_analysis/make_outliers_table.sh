#!/bin/bash

location_of_py_file="make_outliers_table.py"
location_of_data_file="test_input_for_outliers.tsv"
iqrs_over_median=0.5 #Note 1.5 IQRs is suggested, this is just for test data. 
gene_column_name="Genes"
output_file="test_outliers_output_up.tsv"
sample_names_file="test_samples.txt"
aggregate=True
updown="up"

python2.7 ${location_of_py_file} \
--input_df ${location_of_data_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_file ${output_file} \
--sample_names_file ${sample_names_file} \
--aggregate ${aggregate} \
--up_or_down ${updown}
