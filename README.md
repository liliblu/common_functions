# Common_functions
Reusable functions and tools for CPTAC data.  
Some command line tools for easy looping of analysis and use by others.  
Also some sbatch scripts for running jobs on slurm.   

## Outliers analysis  
##### Dependencies:
pandas  
numpy  
scipy.stats  
argparse  

Example below is also in a jupyter notebook in outliers directory. 
Example code for running outliers_takes_nans.py via command line:

```
location_of_py_file="outliers/outliers_takes_nans.py"
location_of_data_file="genes_vs_samples.tsv"
iqrs_over_median=1.5
gene_column_name="geneSymbol"
output_file="outliers_output.tsv"
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


## Group comparisons using outliers

##### Dependencies:
pandas  
numpy  
matplotlib.patches
matplotlib.pyplot  
seaborn  
scipy.stats  
sys  
argparse  
datetime  

Example below is also in a jupyter notebook in outliers directory. 
Example code for running outliers_takes_nans.py via command line:

```
location_of_py_file="outliers/outlier_comparison_generator.py"
location_of_outliers_file="outliers_output.tsv"
gene_column_name="geneSymbol"
fdr_cut_off=0.05
genes_to_highlight="druggable_gene_list.txt"
blue_or_red="red"

output_prefix="prefix_for_output_files"
group1_label="group_of_interest"
group1_list="group_of_interest_samples.txt"
group2_label="not_in_group_of_interest"
group2_list="not_in_group_of_interest_samples.txt"
group_colors="group_colors_for_labels.tsv"

python2.7  ${location_of_py_file} \
--outliers_table  ${location_of_outliers_file} \
--gene_column_name ${gene_column_name} \
--fdr_cut_off ${fdr_cut_off} \
--output_prefix ${output_prefix} \
--group1_label ${group1_label} \
--group1_list ${group1_list} \
--group2_label ${group2_label} \
--group2_list ${group2_list} \
--genes_to_highlight ${genes_to_highlight} \
--blue_or_red ${blue_or_red} \
--group_colors ${group_colors}
```

##### Argument explanation:
*location_of_py_file:* Path to outlier_comparison_generator.py file.  

*location_of_outliers_file:* Path to outliers output from outliers_takes_nans.py.  

*gene_column_name:* Name of column used for labeling rows, usually gene names. 

*fdr_cut_off:* Threshold for signficance, used for selecting genes that are significantly different between two groups, after multiple testing correction.  

*genes_to_highlight:* Optional. Can input a list of genes that will be marked with "\*" in the output heatmap.  

*blue_or_red:* Color map colors for output heatmap. Default red. Useful for differentiating between up and down outliers. 

*output_prefix:* Prefix used for output files. 

*group1_label:* Label used for group of interest. If inputting colors key for groups, must match color key label.  

*group1_list:* Path to list of samples belonging to primary group of interest. Samples not included in this list or group2_list will be ignored in analysis. 

*group2_label:* Label used for samples **not** in group of interest. If inputting colors key for groups, must match color key label.  

*group2_list:* Path to list of samples **not** in group of interest. Samples not included in this list or group1_list will be ignored in analysis.  

*group_colors:* Tab separated list of group labels with paired hex color for label. 

##### Output explanation:
This tool has 2 outputs. The first output is a heatmap, visualizing genes that are significantly enriched in the group of interest, compared to the second group. The second output is a list of genes/labels for rows that were found to be significantly enriched in the group of interest. 

