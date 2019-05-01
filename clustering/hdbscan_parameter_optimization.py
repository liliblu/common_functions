import hdbscan, matplotlib as mpl, seaborn as sns, matplotlib.pyplot as plt
import pandas as pd, numpy as np
from sklearn import metrics, decomposition
import argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--std', type=int)
    parser.add_argument('--nmembers', type=int)
    parser.add_argument('--file_genesets', type=str, default='positive_controls/ptm.sig.db.all.flanking.human.v1.9.0.gmt')
    parser.add_argument('--file_phospho', type=str, default='brca-prospective/brca_prosp_v2.1_phosphoproteome-ratio-norm-NArm.gct_site-centric-seqwin_localized_n130x27352.gct')
    parser.add_argument('--file_samples', type=str, default='brca-prospective/sample_roster.txt')
    parser.add_argument('--geneSymbol_col', type=str, default='geneSymbol')
    parser.add_argument('--psite_col', type=str, default='variableSites')
    parser.add_argument('--output_prefix', type=str, default='')

    args = parser.parse_args()
    std = args.std
    nmembers = args.nmembers
    file_genesets = args.file_genesets
    file_phospho = args.file_phospho
    file_samples = args.file_samples
    geneSymbol_col = args.geneSymbol_col
    psite_col = args.psite_col
    output_prefix = args.output_prefix

    with open(file_genesets, 'r') as fh:
        genesets = {line.split()[0]:line.split()[1:] for line in fh.readlines()}

    reverse = {i:k for k, v in genesets.items() for i in v}
    with open(file_samples, 'r') as fh:
         samples = [x.strip() for x in fh.readlines() if '.' not in x]
    df = pd.read_csv(file_phospho, sep='\t', skiprows=2, index_col=0).replace(['na', 'NA', 'NAN'], np.nan)
    df[psite_col] = df[psite_col].str.strip()
    df = df.loc[(df[geneSymbol_col].isnull()==False)&
                (df[psite_col].str.contains(' ')==False), :]
    df[samples] = df[samples].astype(float)
    df['label'] = df[geneSymbol_col] + '_' + df[psite_col].str[0:-2]

    df = df.loc[(df[samples].isnull().sum(axis=1)<len(samples)/2), ['label']+samples]


    opt_results = pd.DataFrame()

    try:
        temp = pd.read_csv('%sbrca_pros_phospho_corr_std%s.txt'%(output_prefix, std), sep='\t', index_col=0)
    except:

        cut_off = np.nanpercentile(df[samples].std(axis=1), 100-std)
        temp = df.loc[df[samples].std(axis=1)>cut_off, :]
        temp = temp[samples].transpose().corr().abs()
        temp.to_csv('%sbrca_pros_phospho_corr_std%s.txt'%(output_prefix, std), sep='\t')
    try:
        labels = pd.read_csv('%shdbscan_results_%sstd_%s_nmembers.txt'%(output_prefix, std, nmembers), sep='\t')
    except:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=nmembers)
        clusterer.fit(temp)

        labels = pd.DataFrame()
        labels['hdbscan'] = clusterer.labels_
        labels.index = temp.index
        labels = labels.loc[labels['hdbscan']!=-1, :]

        labels = labels.merge(df[['label']], how='left', left_index=True, right_index=True)
        labels['ptm_set'] = [reverse.get(x, np.nan) for x in labels['label']]
        labels.to_csv('%shdbscan_results_%s_std_%s_nmembers.txt'%(output_prefix, std, nmembers), sep='\t')
    labels = labels.loc[labels['ptm_set'].isnull()==False, :]

    genes_left = len(labels)

    adjRand = metrics.adjusted_rand_score(labels['hdbscan'], labels['ptm_set'])
    adjMut = metrics.adjusted_mutual_info_score(labels['hdbscan'], labels['ptm_set'])

    line = {'std':[std], 'nmembers':[nmembers],
            'adjRand':[adjRand], 'adjMut':[adjMut],
            'genesCommon':[genes_left]}

    opt_results = opt_results.append(pd.DataFrame.from_dict(line), ignore_index=True)

    opt_results.to_csv('%sopt_results_%sstd_%snmembers.txt'%(output_prefix, std, nmembers), sep='\t', index=False)
