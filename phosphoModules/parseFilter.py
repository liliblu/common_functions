import pandas as pd
import numpy as np

true_label_col = 'true_label'
clusterlabels_col = 'labels'

gmt_site_sep = '_'
site_id_sep = '-'
site_id_col = 'site_id'
psite_sep = ' '
numNA_col = 'numNA'
frac_noNA = 0.5

def parseFilterGCT(gct_file, gene_col, psite_col, samples):
    if gct_file[-4:] == '.gct':
        df = pd.read_csv(gct_file, sep='\t', skiprows=2, low_memory=False)
    else:
        df = pd.read_csv(gct_file, sep='\t', low_memory=False)
    df = df.replace(['na', 'NaN', 'nan', 'NA', 'NAN', 'Nan'],
                    np.nan).dropna(subset=[gene_col], axis=0).astype(float, errors='ignore')
    df[site_id_col] = df[gene_col] + site_id_sep + df[psite_col]
    df = df.loc[df[samples].isnull().sum(axis=1) < len(samples) * frac_noNA, :]
    return df

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
                    subset = subset.groupby(site_id_col).agg(lambda row: row.mean())
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)

                elif len(subset) == 1:
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)


            elif len(subset) == 1:
                deduped_phospho = deduped_phospho.append(subset, ignore_index=True)
    deduped_phospho[site_id_col] = deduped_phospho[site_id_col].str.replace(r'[a-z]', '')
    deduped_phospho = deduped_phospho.set_index(site_id_col)
    return deduped_phospho

def convertLineToResiduals(ph, prot):
    nonull = ((ph.isnull() == False) & (prot.loc[ph.name[0], :].isnull() == False))

    features = prot.loc[ph.name[0], :][nonull].values.reshape(-1, 1)
    labels = ph[nonull].values
    model = lm.LinearRegression().fit(features, labels)
    prediction = model.predict(features)
    residuals = labels - prediction
    return pd.Series(residuals, index=nonull[nonull].index)

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
    with open(gmt_file, 'r') as fh:
        genesets = {line.split()[0]: [site.split(':')[0] for site in line.split()[1].split('|')[1:]] for line in
                    fh.readlines()}
    return genesets

def filterSTDcorr(df, samples, std):
    cut_off = np.nanpercentile(df[samples].std(axis=1), 100 - std)
    df = df.loc[((df[samples].std(axis=1))>cut_off), samples]
    corr = df.transpose().corr().abs()
    return corr
