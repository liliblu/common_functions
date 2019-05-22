import hdbscan
import pandas as pd
import numpy as np
from sklearn import metrics

true_label_col = 'true_label'
clusterlabels_col = 'labels'
site_id_sep = '-'
site_id_col = 'site_id'
psite_sep = ' '
numNA_col = 'numNA'
frac_noNA = 0.5


def collectSites(df):
    gene_sites = {}
    for label in df[site_id_col]:
        gene = label.split(site_id_sep)[0]
        sites = label.split(site_id_sep)[1].strip().split(psite_sep)
        gene_sites[gene] = gene_sites.get(gene, []) + sites
    gene_sites = {gene:set(sites) for gene, sites in gene_sites.items()}
    return gene_sites


def removeDups(df, gene_col, psite_col, samples):
    gene_sites = collectSites(df)
    cols = [gene_col, psite_col, site_id_col, numNA_col]
    # sites_picked_by_na = []
    # sites_picked_randomly = []
    # sites_unique = []
    df[numNA_col] = df[samples].isnull().sum(axis=1)
    deduped_phospho = pd.DataFrame()
    for gene, sites in gene_sites.items():
        for site in sites:
            gene_site = gene + site_id_sep + site
            subset = df.loc[((df[gene_col] == gene) &
                             (df[psite_col].str.contains(site))), cols + samples]
            subset[site_id_col] = gene_site

            if len(subset) > 1:
                subset = subset.loc[subset[numNA_col] == subset[numNA_col].min(), :]
                if len(subset) > 1:
                    subset = subset.sample(1)
                    # sites_picked_randomly.append(gene_site)
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)

                elif len(subset) == 1:
                    # sites_picked_by_na.append(gene_site)
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)


            elif len(subset) == 1:
                # sites_unique.append(gene_site)
                deduped_phospho = deduped_phospho.append(subset, ignore_index=True)
    deduped_phospho[site_id_col] = deduped_phospho[site_id_col].str.replace(r'[a-z]', '')
    deduped_phospho = deduped_phospho.set_index(site_id_col)
    return deduped_phospho


def parseFilterGCT(gct_file, gene_col, psite_col, samples):
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, low_memory=False)
    df = df.replace(['na', 'NaN', 'nan', 'NA', 'NAN', 'Nan'],
                    np.nan).dropna(subset=['geneSymbol'], axis=0).astype(float, errors='ignore')
    df[site_id_col] = df[gene_col] + site_id_sep + df[psite_col]
    df = df.loc[df[samples].isnull().sum(axis=1) < len(samples) * frac_noNA, :]
    return df


def preprocessGCT(gct_file,
                  samples,
                  gene_col = None,
                  psite_col = None,
                  writeToFilePrefix = None):
    if gene_col is None:
        gene_col = 'geneSymbol'
    if psite_col is None:
        psite_col = 'variableSites'
    try:
        df = pd.read_csv(writeToFilePrefix + '.preprocessed.txt', sep='\t', index_col=0)
    except:
        df = parseFilterGCT(gct_file, gene_col, psite_col, samples)
        df = removeDups(df, gene_col, psite_col, samples)

    if not (writeToFilePrefix is None):
        df.to_csv(writeToFilePrefix + '.preprocessed.txt', sep='\t')
    return df


def parseGMT(gmt_file):
    genesets = {}
    with open(gmt_file, 'r') as fh:
        genesets = {line.split()[0]: [site.split(':')[0] for site in line.split()[1].split('|')[1:]] for line in
                    fh.readlines()}
    return genesets


def evalutateLabelsRand(labelsDF, genesets):
    score = 0
    sets_compared = 0
    for setName, sites in genesets.items():
        sites = [site.replace('_', '-') for site in sites]
        if any([site in labelsDF.index for site in sites]):
            sets_compared += 1
            labelsDF[true_label_col] = labelsDF.index.isin(sites)
            score += metrics.adjusted_rand_score(labelsDF[true_label_col], labelsDF[clusterlabels_col])
    return score, sets_compared


def evalulateLabelsSilhouette(df, labelDF):
    df = df.reindex(labelDF.index, axis=0).reindex(labelDF.index, axis=1)
    score = metrics.silhouette_score(df, labelDF[clusterlabels_col])
    return score


def filterSTDcorr(df, samples, std):
    cut_off = np.nanpercentile(df[samples].std(axis=1), 100 - std)
    df = df.loc[((df[samples].std(axis=1))>cut_off), samples]
    corr = df.transpose().corr().abs()
    return corr


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
                         nmembers_range: list = None) -> dict:
    positive_controls = not (gmt_file is None)
    if std_range is None:
        std_range = range(5, 40, 5)
    if nmembers_range is None:
        nmembers_range = range(2, 11)

    scores = pd.DataFrame()
    for std in std_range:
        corr = filterSTDcorr(df, samples, std).dropna(how='any', axis=0)
        corr = corr.reindex(corr.index, axis=1)
        for nmembers in nmembers_range:
            line = {"nmembers":[nmembers], "std":[std], "nSites":[len(corr)]}
            labels = runHdbscan(corr, nmembers)
            if positive_controls:
                genesets = parseGMT(gmt_file)
                rand, sets_compared = evalutateLabelsRand(labels, genesets)
                line.update({"adjRand":[rand], "randSetsCompared":[sets_compared]})

            silhouette = evalulateLabelsSilhouette(corr, labels)
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
             nmembers_range=None):
    if std_range is None:
        std_range = range(5, 40, 5)
    if nmembers_range is None:
        nmembers_range = range(2, 11)

    with open(samples_file, 'r') as fh:
        samples = [x.strip() for x in fh.readlines()]

    df = preprocessGCT(gct_file, samples,
                       gene_col=gene_col,
                       psite_col=psite_col,
                       writeToFilePrefix=writeToFilePrefix)

    scores = runThroughParameters(df, samples,
                                  gmt_file=gmt_file,
                                  std_range=std_range,
                                  nmembers_range=nmembers_range)

    return scores
