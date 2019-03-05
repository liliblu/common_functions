#!/bin/bash

location_of_py_file="../../../common_functions/adjusted_outliers.py"
location_of_outliers_file="../data/UCEC_phosphoproteomics_pass_only.cct"
iqrs_over_median=1.5
gene_column_name="geneSymbol"
output_file="test.txt"
sample_names_file="../samples_endo.txt"
updown="up"

python2.7 ${location_of_py_file} \
--input_df ${location_of_outliers_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_file ${output_file} \
--sample_names_file ${sample_names_file} \
--up_or_down ${updown}
