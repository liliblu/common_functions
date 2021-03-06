from typing import Optional, Iterable
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scipy.stats
import catheat
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
import sklearn.linear_model as lm
import warnings

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
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
        for i in range(0, int(n)-1):
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

def listToFile(lis, file_name):
    with open(file_name, 'w') as fh:
        for x in lis:
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
        row_order = df.index

    if col_cluster==True:
        col_order = computeOrder(df.transpose(), optimal, dist_method, cluster_method)
        col_order = [df.columns[i] for i in col_order]
    else:
        col_order = df.columns

    df = df.reindex(col_order, axis=1).reindex(row_order, axis=0)

    ax = sns.heatmap(df, **heatmap_kws)

    return ax, df

def convertLineToResiduals(ph, prot, alphas=[2**i for i in range(-10, 10, 1)], cv=5):
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    nonull = ((ph.isnull() == False) & (prot.loc[ph.name[0], :].isnull() == False))
    if sum(nonull) < cv:
        residuals = np.empty(len(ph))
        residuals = pd.Series(residuals, index=ph.index)
    else:
        features = prot.loc[ph.name[0], :][nonull].values.reshape(-1, 1)
        labels = ph[nonull].values
        model = lm.RidgeCV(alphas=alphas, cv=cv).fit(features, labels)
        prediction = model.predict(features)
        residuals = labels - prediction
        residuals = pd.Series(residuals, index=ph[nonull].index)
    return residuals


def corrNA(array1, array2):
    if isinstance(array1, pd.Series):
        nonull = ((array1.isna()==False) & (array2.values.isna()==False))
        return scipy.stats.pearsonr(array1[nonull], array2[nonull])
    nonull = ((np.isnan(array1)==False) & (np.isnan(array2)==False))
    return scipy.stats.pearsonr(array1[nonull], array2[nonull])


def MAFsampleName(string):
    string = string.split('_')[0]
    if 'CPT' in string:
        return string
    else:
        return 'X'+string


def deduplicate_rows(
    df: pd.DataFrame,
    samples: Iterable,
    sequence_col: Optional[str] = 'sequence',
    maximize_cols: Iterable = ['bestScore', 'Best_scoreVML', 'bestDeltaForwardReverseScore']
    ) -> pd.DataFrame:

    groupby = df.index.names
    df = df.reset_index()

    df[maximize_cols] = df[maximize_cols].astype(float)
    df['numNAs'] = df[samples].isnull().sum(axis=1)
    if sequence_col:
        df['len'] = df[sequence_col].str.len()
        return df.groupby(groupby).apply(lambda row: row.nsmallest(1, columns=['numNAs','len'], keep='all').nlargest(1, columns=maximize_cols, keep='first')).reset_index(level=-1, drop=True)
    return df.groupby(groupby).apply(lambda row: row.nsmallest(1, columns=['numNAs'], keep='all').nlargest(1, columns=maximize_cols, keep='first')).reset_index(level=-1, drop=True)


def parseGCT(
    filename: str,
    typ: str = 'phospho',
    inds: Iterable = ['geneSymbol'],
    sampleid_row_name: str = 'Sample.ID',
    writefile: bool = False,
    output_prefix: Optional[str] = "deduped",
    **dedup_kws
    ) -> pd.DataFrame:

    df = pd.read_csv(filename, sep='\t', skiprows=2, low_memory=False).replace(['na', 'NA', 'Na', 'nan', 'NAN', 'NaN', 'Nan'], np.nan)
    samples = list(df.loc[df[df.columns[0]]==sampleid_row_name, :].dropna(axis=1))[1:]

    df = df.dropna(subset=inds, how='all')
    df = df.set_index(inds)
    df[samples] = df[samples].astype(float, errors='ignore')

    if df.index.duplicated().sum() > 0:
        df = deduplicate_rows(df, samples, **dedup_kws)
    df = df[samples]

    if writefile:
        df.to_csv('%s.csv'%output_prefix)
    return df
