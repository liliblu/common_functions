import pandas as pd
import numpy as np
from sklearn import metrics

site_id_sep = '-'

def evalutateLabelsRand(labelsDF, genesets):
    score = 0
    sets_compared = 0
    for setName, sites in genesets.items():
        sites = [site.replace('_', site_id_sep) for site in sites]
        if any([site in labelsDF.index for site in sites]):
            sets_compared += 1
            labelsDF[true_label_col] = labelsDF.index.isin(sites)
            score += metrics.adjusted_rand_score(labelsDF[true_label_col], labelsDF[clusterlabels_col])
    return score, sets_compared


def evalulateLabelsSilhouette(df, labelDF):
    df = df.reindex(labelDF.index, axis=0).reindex(labelDF.index, axis=1)
    score = metrics.silhouette_score(df, labelDF[clusterlabels_col])
    return score
