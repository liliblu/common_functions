from __future__ import division

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scipy.stats
import catheat

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


def makeRankedHeatmaps(input_df, sampsIn, sampsOut, prefix, colors):
    #Useful for showing RRA analysis
    num_plots = len(input_df)
    num_samples = len(list(input_df))
    if num_plots > 0:
        fig, axs = plt.subplots(nrows=num_plots, ncols=1,
                                figsize=(0.08*num_samples, 0.2*num_plots),
                                gridspec_kw=dict(hspace=0.05))

        ax_ind = 0
        for i, row in input_df.iterrows():
            row = row.sort_values(na_position='first')
            nas = row.isnull()
            ax=axs[ax_ind]

            row[sampsIn] = 'In group'
            row[sampsOut] = 'Out of group'
#             row[unknown_samps] = 'Unknown'
            row[nas] = 'Na'

            row = pd.DataFrame(row).transpose()
            catheat.heatmap(row, cmap=colors,
                            xticklabels=False,
                            yticklabels=False,
                            ax=ax,
                            legend=False,
                           )
            ax.set_ylabel(i, rotation=0, fontsize=10, ha="right", va='center')

            ax_ind+=1

        axs[0].set_title('Phosphosites ranked in %s' %prefix, fontsize=16)
        plt.text(0, 1.3, 'Low', fontsize=16, ha='left', va='bottom', transform=axs[0].transAxes)
        plt.text(1, 1.3, 'High', fontsize=16, ha='right', va='bottom', transform=axs[0].transAxes)

        handles = [mpatches.Patch(color=color, label=label) for label, color in colors.iteritems()]
        plt.legend(handles=handles, loc=(1.01, 0))

    return axs

def listToFile(list, file_name):
    with open(file_name, 'w') as fh:
        for x in list:
            fh.write('%s\n'%x)

def computeOrder(df,
                 optimal=True,
                 dist_method="euclidean",
                 cluster_method="average"):

    dist_mat = pdist(df, metric=dist_method)
    link_mat = hierarchy.linkage(dist_mat, method=cluster_method)

    if optimal==True:
        return hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(link_mat, dist_mat))
    else:
        return hierarchy.leaves_list(link_mat)

def clustermap(df,
               dist_method="euclidean",
               cluster_method="average",
               col_cluster=True,
               row_cluster=True,
               optimal=True,
               **heatmap_kws):

    if row_cluster==True:
        row_order = computeOrder(df, optimal, dist_method, cluster_method)
        row_order = [df.index[i] for i in row_order]
    else:
        col_order = df.index

    if col_cluster==True:
        col_order = computeOrder(df.transpose(), optimal, dist_method, cluster_method)
        col_order = [df.columns[i] for i in col_order]
    else:
        row_order = df.columns

    df = df.reindex(col_order, axis=1).reindex(row_order, axis=0)

    ax = sns.heatmap(df, **heatmap_kws)

    return ax, df
    
