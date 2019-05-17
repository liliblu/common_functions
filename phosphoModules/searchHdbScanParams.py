import hdbscan
import pandas as pd
import numpy as np
from sklearn import metrics

true_label_col = 'true_label'
clusterlabels_col = 'labels'

def parseGMT(gmt_file):
    genesets = {}
    with open(gmt_file, 'r') as fh:
        genesets = {line.split()[0]:[site.split(':')[0] for site in line.split()[1].split('|')[1:]] for line in fh.readlines()}
    return genesets

def evalutateLabelsRand(labelsDF, genesets, site_id_col='site_id'):
    score = 0
    for setName, sites in genesets.items():
        labelsDF[true_label_col] = labelsDF[site_id_col].isin(sites)
        score+=metrics.adjusted_rand_score(labelsDF[true_label_col], labelsDF[clusterlabels_col])
    return score

def evalulateLabelsSilhouette(df, labelDF):
    df = df.reindex(labelDF.index, axis=0).reindex(samples, axis=1)
    score = metrics.silhouette_score(df, labelDF[clusterlabels_col].values())
    return score

def filterSTDcorr(df, samples, std):
    cut_off = np.nanpercentile(df[samples].std(), 100-std)
    df = df.loc[df[samples].std()>cut_off, samples]
    corr = df.transpose().corr().abs()
    return corr

def runHdbscan(corr, nmembers):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=nmembers)
    clusterer.fit(corr)
    labels = pd.DataFrame()
    labels[clusterlabels_col] = clusterer.labels_
    labels.index = corr.index
    labels = labels.loc[labels[clusterlabels_col]!=-1, :]
    return labels


def runThroughParameters(df, samples, gmt_file=None,
                         positive_controls=False,
                         std_range=range(5, 40, 5),
                         nmembers_range=range(2, 11)):
    scores = {'adjRand':{}, 'silhouette':{}}
    for std in std_range:
        corr = filterSTD(df, samples, std)
        for nmembers in nmembers_range:
            key = (std, nmembers)
            labels = runHdbscan(corr, nmembers)
            if positive_controls:
                genesets = parseGMT(gmt_file)
                scores['adjRand'][key] = evalutateLabelsRand(labels, genesets)
            
            scores['silhouette'][key] = evalulateLabelsSilhouette(corr, labels)

    return scores
