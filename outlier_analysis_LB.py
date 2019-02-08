from __future__ import division
import pandas as pd
import numpy as np
import scipy.stats
import argparse

# to do: convert na to Nan, convert to float.
def cleanDF(df, sample_columns):
    '''
    Convert string nans to np.nan and string numbers to floats.
    '''
    df = df.replace(['na', 'NaN', 'Na', 'nan'], np.nan)
    df[sample_columns] = df[sample_columns].astype(float)

    return df

def getColumns(df, sample_columns):
    '''
    Pulls out annotation columns. For phospho data there should be 3,
    (gene symbol, isoform ID, phosphosite), for proteome data there should be 2,
    (gene symbol, isoform ID).
    Also generates a list of outlier names, needed for convertToOutliers.
    '''
    gene_info_columns = list(set(list(df)) - set(sample_columns))

    return gene_info_columns


def convertToOutliers(df, gene_info_columns, sample_columns, NUM_IQRs, up_or_down):
    '''
    Calculates the median, and inter-quartile range for each row/isoform.
    Inter-quartile range is defined as the value difference between the 25th and 75th percentile.
    Here, NaNs are ignored for each row, therefore a different number of values may be used for each row.

    '''
    df['row_iqr'] = scipy.stats.iqr(df[sample_columns], axis=1, nan_policy='omit')
    df['row_median'] = np.nanquantile(df[sample_columns], q=0.5, axis=1)

    df['row_medPlus'] = (df['row_median'] + (NUM_IQRs*df['row_iqr']))
    df['row_medMinus'] = (df['row_median'] - (NUM_IQRs*df['row_iqr']))

    outlier_df = pd.DataFrame()
    outlier_df[gene_info_columns] = df[gene_info_columns]
    outlier_df = outlier_df.assign(**{col:0 for col in sample_columns})

    if up_or_down == 'up':
        for col in sample_columns:
            outlier_df.loc[(df[col] > df['row_medPlus']), col] = 1
            outlier_df.loc[(df[col].isnull()), col] = np.nan

    elif up_or_down == 'down':
        for col in sample_columns:
            outlier_df.loc[(df[col] > df['row_medPlus']), col] = 1
            outlier_df.loc[(df[col].isnull()), col] = np.nan

    elif up_or_down == 'both':
        for col in sample_columns:
            outlier_df.loc[(df[col] > df['row_medPlus']), col] = 1
            outlier_df.loc[(df[col] < df['row_medMinus']), col] = -1
            outlier_df.loc[(df[col].isnull()), col] = np.nan

    return outlier_df


def aggregatePhosphosites(df, isoform_column_name, sample_columns):
    df.loc[:,'num_psites'] = np.nan
    count_df = df.groupby(by=isoform_column_name)['num_psites'].agg(lambda x: len(x)).reset_index(name='counts')
    df = df.groupby(by=isoform_column_name)[sample_columns].agg((lambda x: sum(x)))
    df = df.merge(count_df, how='left', on=isoform_column_name)

    return df

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--input_df', type=str)
    parser.add_argument('--experiment_type', type=str, choices=['phospho', 'not_phospho'], default='not_phospho')
    parser.add_argument('--iqrs_over_median', type=float, default=1.5)
    parser.add_argument('--isoform_column_name', type=str, default='id')
    parser.add_argument('--output_file', type=str, default='outliers.tsv')
    parser.add_argument('--sample_names_file', type=str, default='sample_roster.txt')
    parser.add_argument('--up_or_down', type=str, choices=['up', 'down', 'both'], default='up')

    args = parser.parse_args()

    data_input = args.input_df
    experiment_type = args.experiment_type
    isoform_column_name = args.isoform_column_name
    write_results_to = args.output_file
    NUM_IQRs = args.iqrs_over_median
    sample_names = args.sample_names_file
    up_or_down = args.up_or_down

    sample_columns = []
    with open(sample_names, 'r') as names_file:
        for name in names_file.readlines():
            name = name.strip()
            sample_columns.append(name)

    sample_data = pd.read_csv(data_input, sep='\t')

    sample_data = cleanDF(sample_data, sample_columns)

    gene_info_columns = getColumns(sample_data, sample_columns)

    outliers = convertToOutliers(sample_data, gene_info_columns, sample_columns, NUM_IQRs, up_or_down)

    if experiment_type == 'phospho':
        outliers = aggregatePhosphosites(outliers, isoform_column_name, sample_columns)
        outliers.to_csv(write_results_to, sep='\t', index=False)
    elif experiment_type == 'not_phospho':
        outliers[[isoform_column_name] + sample_columns].to_csv(write_results_to, sep='\t', index=False)

    print('Outlier analysis complete. Results are in %s' %write_results_to)
