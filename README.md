# Common_functions
Reusable functions and tools for CPTAC data.  
Some command line tools for easy looping of analysis and use by others.  
Also some sbatch scripts for running jobs on slurm.   

##### Dependencies:
pandas
numpy
scipy.stats
argparse

Below examples are also in a jupyter notebook in outliers directory. 
Example code for running outliers_takes_nans.py via command line:

```
location_of_py_file="common_functions/adjusted_outliers.py"
location_of_data_file="genes_vs_samples.tsv"
iqrs_over_median=1.5
gene_column_name="geneSymbol"
output_file="outliers.tsv"
sample_names_file="columns_to_consider.txt"
updown="up"

python2.7 ${location_of_py_file} \
--input_df ${location_of_data_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_file ${output_file} \
--sample_names_file ${sample_names_file} \
--up_or_down ${updown}
```

##### Argument explanation:
*location_of_py_file:* Path to outliers_takes_nan.py file. 

*location_of_data_file:* Path to input file. Sample names must match column names. Needs gene names column as well. Can also include other annotation columns, but those will not be propagated to output.  

*iqrs_over_median:* Threshold for a value to be considered an "outlier". Number of inter-quartile ranges over the median, calculated for each gene, to be used as a threshold.  

*gene_column_name:* The name of the column to be used for aggregating values per gene. For instance, if phospho-peptide data is used as input, it may be helpful to aggregate outliers at the gene or isoform level.  

*output_file:* Name/location to put outliers output.   

*sample_names_file:* Column/sample names to be used for analysis. Columns not included in this list will be ignored.   

*up_or_down*: Whether analysis should mark outliers that are **greater than** median + IQR threshold, or **less than** median - IQR threshold.  

##### Output explanation:
The output will have a column corresponding to the gene_column_name from the input. For each sample listed in sample_names_file, there will be 2 output columns. For sample X, one column will be called X_outliers which counts the number of outliers for that gene/row for that sample. The second columns will be called X_notOutliers which counts the number of sites per gene/row that are not outliers in that sample. For input data where each each value in gene_column_name is unique, and doesn't have missing values, the values from X_outliers and X_notOutliers will always add up to 1. In other words, each sample will be either an outlier or not an outlier (or be a missing value). This format is helpful for calculating comparisons in downstream analysis. 
