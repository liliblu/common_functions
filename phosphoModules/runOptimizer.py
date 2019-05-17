from preprocessing import *
from searchHdbScanParams import *
import hdbscan
import pandas as pd
import numpy as np
from sklearn import metrics

# site_id_sep = '-'
# site_id_col = 'site_id'
# psite_sep = ' '
# numNA_col = 'numNA'
# true_label_col = 'true_label'
# clusterlabels_col = 'labels'

def optimize(gct_file,
             samples_file,
             gmt_file=None,
             gene_col='geneSymbol',
             psite_col='variableSites',
             writeToFilePrefix=None,
             std_range=range(5, 40, 5),
             nmembers_range=range(2, 11)):

    with open(samples_file, 'r') as fh:
        samples = [x.strip() for x in fh.readlines()]

    df = preprocessGCT(gct_file, samples,
                       gene_col=gene_col,
                       psite_col=psite_col,
                       writeToFilePrefix=writeToFilePrefix)

    positive_controls = gmt_file!=None
    scores = runThroughParameters(df, samples,
                                gmt_file=gmt_file,
                                positive_controls=positive_controls,
                                std_range=std_range,
                                nmembers_range=nmembers_range)
    return scores
