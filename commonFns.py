from __future__ import division

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scipy.stats


def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = sum(~np.isnan(pvalues))
    new_pvalues = empty(len(pvalues))

    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues

    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values = [x for x in values if ~np.isnan(x[0])]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":

        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values = [x for x in values if ~np.isnan(x[0])]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]

    new_pvalues[np.isnan(pvalues)] = np.nan
    for i, val in enumerate(new_pvalues):
        if ~np.isnan(val):
            new_pvalues[i] = min(1, val)

    return new_pvalues

def makeHeatMapTable(pval_table, outlier_table, pval_column, sig_threshold, samples_in_subtype, samples_out_subtype, isoform_column='id'):
    sig_isoforms = list(pval_table.loc[(pval_table[pval_column] <= sig_threshold), isoform_column])
    heatmap_table = outlier_table.loc[(outlier_table[isoform_column].isin(sig_isoforms) == True), [isoform_column] + samples_in_subtype + samples_out_subtype]
    heatmap_table = heatmap_table.set_index(isoform_column)
    return heatmap_table


def fileToList(group_list):
    with open(group_list, 'r') as fh:
        return [line.strip() for line in fh.readlines()]

def fileToDict(tsv_map_file_name):
    with open(tsv_map_file_name, 'r') as fh:
        return {line.split()[0]:line.split()[1] for line in fh.readlines()}

def compareQuantiles(row, order, fraction):
    quartile_length = int(len(order)*fraction)
    top_quartile = row[order[0:quartile_length]]
    bottom_quartile = row[order[-quartile_length:]]

    top_outliers = top_quartile.sum()
    top_not_outliers = (row['counts']*quartile_length*2) - top_outliers

    bottom_outliers = bottom_quartile.sum()
    bottom_not_outliers = (row['counts']*quartile_length*2) - bottom_outliers

    fisher_pval = scipy.stats.fisher_exact(np.array([[top_outliers, top_not_outliers], [bottom_outliers, bottom_not_outliers]]))[1]

    return fisher_pval

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def geometric(x, *p):
    A, pr = p
    return A*np.exp(((1-pr)**(x-1))*pr)

def gaussPlusGeo(x, *p):
    gA, gmu, gsigma, geoA, geopr = p
    gaussp = [gA, gmu, gsigma]
    geop = [geoA, geopr]
    return gauss(x, *gaussp) + geometric(x, *geop)
