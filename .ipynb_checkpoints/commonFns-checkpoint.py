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
    n = len(pvalues)
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
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
    return new_pvalues


def testDifferentGroupsOutliers(sample_set1, sample_set2, outlier_table, psite_count_column='counts', phospho=True):

    if phospho == True:
        outlier_table['Outlier1'] = outlier_table[sample_set1].sum(axis=1)
        outlier_table['NotOutlier1'] = ((outlier_table[psite_count_column]*len(sample_set1)) - outlier_table['Outlier1'])

        outlier_table['Outlier2'] = outlier_table[sample_set2].sum(axis=1)
        outlier_table['NotOutlier2'] = ((outlier_table[psite_count_column]*len(sample_set2)) - outlier_table['Outlier2'])
    elif phospho == False:
        outlier_table['Outlier1'] = (outlier_table[sample_set1] > 0).sum(axis=1)
        outlier_table['NotOutlier1'] = (outlier_table[sample_set1] == 0).sum(axis=1)

        outlier_table['Outlier2'] = (outlier_table[sample_set2] > 0).sum(axis=1)
        outlier_table['NotOutlier2'] = (outlier_table[sample_set2] == 0).sum(axis=1)

    outlier_table['fisherp'] = outlier_table.apply((lambda r: scipy.stats.fisher_exact([[r['Outlier1'],
                                                                                     r['Outlier2']],
                                                                                    [r['NotOutlier1'],
                                                                                     r['NotOutlier2']]])[1]),axis=1)

    outlier_table['fisherFDR'] = correct_pvalues_for_multiple_testing(list(outlier_table['fisherp']),
                                     correction_type = "Benjamini-Hochberg")

    return outlier_table['fisherFDR']


def makeHeatMapTable(pval_table, outlier_table, pval_column, sig_threshold, samples_in_subtype, samples_out_subtype, isoform_column='id'):
    sig_isoforms = list(pval_table.loc[(pval_table[pval_column] <= sig_threshold), isoform_column])
    heatmap_table = outlier_table.loc[(outlier_table[isoform_column].isin(sig_isoforms) == True), [isoform_column] + samples_in_subtype + samples_out_subtype]
    heatmap_table = heatmap_table.set_index(isoform_column)
    return heatmap_table


def makeHeatMap(heatmap_table,
                samples_in_subtype, samples_out_subtype,
                subtype_colors,
                isoform_annotations,
                subtype,
                alt,
                fig_name,
                fig_title = 'Differentially phosphorylated isoforms',
                isoform_column='id',
                yticklabels=False,
                dendrogram=False
                ):

    lut = {sample:subtype_colors[subtype] for sample in samples_in_subtype}
    lut2 = {sample:subtype_colors[alt] for sample in samples_out_subtype}
    lut.update(lut2)

    group = heatmap_table.columns
    column_colors = group.map(lut)

    if yticklabels == True:
        labels = [isoform_annotations.loc[(isoform_annotations[isoform_column] == x), 'geneSymbol'].values[0] for x in heatmap_table.index]
    elif yticklabels == False:
        labels = False

    g = sns.clustermap(heatmap_table,
#                        cmap = sns.light_palette('red'),
                       cmap=sns.cubehelix_palette(start=0.857, rot=0.00, gamma=1.5, hue=1, light=1, dark=0.2, reverse=False, as_cmap=True),
                       col_cluster = False,
                       col_colors = column_colors,
                       xticklabels=False,
                       yticklabels=labels,
#                        z_score = 0,
                       cbar_kws={'label':'# outliers',
#                                  'ticks':np.arange(0, heatmap_table.values.max())
                                },
                       vmin=0,
                       vmax=np.percentile(heatmap_table.values.max(), 70)
                      )
    g.ax_row_dendrogram.set_visible(dendrogram)
    g.cax.set_position((0.15,0.12,0.04,0.6)) #move colorbar to right
    ax = g.ax_heatmap
    ax.set_ylabel('') #change the gene label

#     this chunk makes the legend the describes the different sample groups
    subtype_patch = mpatches.Patch(color=subtype_colors[subtype], label=subtype)
    other_patch = mpatches.Patch(color=subtype_colors['Other'], label='Other')
    fig = plt.gcf()
    fig.legend(handles=[subtype_patch, other_patch], bbox_to_anchor=(0.6, 0.11))

#     give the plot a title
    fig.text(0.4, 0.75, fig_title, fontsize=16)

    #save the plot
    plt.savefig(fig_name, dpi=500, bbox_inches='tight', pad_inches=1)
    plt.show()
    plt.close()


def fileToList(group_list):
    group = []
    with open(group_list, 'r') as fh:
        for line in fh.readlines():
            group.append(line.strip())

    return group

def fileToDict(tsv_map_file_name):
    map = dict()
    with open(tsv_map_file_name, 'r') as fh:
        for line in fh.readlines():
            map[line.split()[0]] = line.split()[1]
    return map


def compareOutliers(row, order, fraction):
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