import hdbscan
import pandas as pd
import parseFilter
import metrics

true_label_col = 'true_label'
clusterlabels_col = 'labels'

gmt_site_sep = '_'
site_id_sep = '-'
site_id_col = 'site_id'
psite_sep = ' '
numNA_col = 'numNA'
frac_noNA = 0.5

def runHdbscan(corr, nmembers):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=nmembers)
    clusterer.fit(corr)
    labels = pd.DataFrame()
    labels[clusterlabels_col] = clusterer.labels_
    labels.index = corr.index
    labels = labels.loc[labels[clusterlabels_col] != -1, :]
    return labels


def runThroughParameters(df, samples,
                         gmt_file: str = None,
                         std_range: list = None,
                         nmembers_range: list = None,
                         write_labels:bool = False,
                         writeToFilePrefix: str = None,
                         ) -> dict:
    positive_controls = not (gmt_file is None)
    if std_range is None:
        std_range = range(5, 40, 5)
    if nmembers_range is None:
        nmembers_range = range(2, 11)

    scores = pd.DataFrame()
    for std in std_range:
        corr = parseFilter.filterSTDcorr(df, samples, std).dropna(how='any', axis=0)
        corr = corr.reindex(corr.index, axis=1)
        for nmembers in nmembers_range:
            labels = runHdbscan(corr, nmembers)
            if write_labels:
                labels.to_csv(writeToFilePrefix+'_labels_%sstd_%snmem.txt'%(std, nmembers), sep='\t')
            line = {"nmembers":[nmembers], "std":[std], "nSites":[len(labels)]}
            if positive_controls:
                genesets = parseGMT(gmt_file)
                rand, sets_compared = metrics.evalutateLabelsRand(labels, genesets)
                line.update({"adjRand":[rand], "randSetsCompared":[sets_compared]})

            silhouette = metrics.evalulateLabelsSilhouette(corr, labels)
            line.update({"sillhouette":[silhouette]})
            scores = scores.append(pd.DataFrame.from_dict(line), ignore_index=True)

    return scores


def optimize(gct_file,
             samples_file,
             gmt_file=None,
             gene_col=None,
             psite_col=None,
             writeToFilePrefix=None,
             std_range=None,
             nmembers_range=None,
             write_labels=False):
    if std_range is None:
        std_range = range(5, 40, 5)
    if nmembers_range is None:
        nmembers_range = range(2, 11)

    with open(samples_file, 'r') as fh:
        samples = [x.strip() for x in fh.readlines()]

    df = parseFilter.preprocessGCT(gct_file, samples,
                       gene_col=gene_col,
                       psite_col=psite_col,
                       writeToFilePrefix=writeToFilePrefix)

    scores = runThroughParameters(df, samples,
                                  gmt_file=gmt_file,
                                  std_range=std_range,
                                  nmembers_range=nmembers_range,
                                  write_labels=write_labels,
                                  writeToFilePrefix=writeToFilePrefix)

    return scores
