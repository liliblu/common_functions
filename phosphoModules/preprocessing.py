import pandas as pd
import numpy as np

site_id_sep = '-'
site_id_col = 'site_id'
psite_sep = ' '
numNA_col = 'numNA'

def collectSites(df):
    gene_sites = {}
    for label in df[site_id_col]:
        gene = label.split(site_id_sep)[0]
        sites = label.split(site_id_sep)[1].strip().split(psite_sep)
        gene_sites[gene] = gene_sites.get(gene, [])+sites
    return gene_sites

def removeDups(df, gene_col, psite_col, samples):
    gene_sites = collectSites(df)
    cols = [gene_col, psite_col, site_id_col]
    sites_picked_by_na = []
    sites_picked_randomly = []
    sites_unique = []
    df[numNA_col] = df[samples].isnull().sum(axis=1)
    for gene, sites in gene_sites.items():
        for site in sites:
            gene_site = gene+site_id_sep+site
            subset = df.loc[((df[gene_col]==gene)&
                             (df[psite_col].str.contains(site))), cols+samples]
            subset[site_id_col] = gene_site

            if len(subset)>1:
                subset = subset.loc[subset[numNA_col]==subset[numNA_col].min(), :]
                if len(subset)>1:
                    subset = subset.sample(1)
                    sites_picked_randomly.append(gene_site)
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)

                elif len(subset)==1:
                    sites_picked_by_na.append(gene_site)
                    deduped_phospho = deduped_phospho.append(subset, ignore_index=True)


            elif len(subset)==1:
                sites_unique.append(gene_site)
                deduped_phospho = deduped_phospho.append(subset, ignore_index=True)
    deduped_phospho[site_id_col] = deduped_phospho[site_id_col].str.replace(r'[a-z]', '')
    deduped_phospho = deduped_phospho.set_index(site_id_col)
    return deduped_phospho

def parseGCT(gct_file, gene_col, psite_col):
    df = pd.read_csv(gct_file, sep='\t', skiprows=2)
    df = df.replace(['na', 'NaN', 'nan', 'NA', 'NAN', 'Nan'],
                    np.nan).dropna(subset=['geneSymbol'], axis=0).astype(float, errors='ignore')
    df[site_id_col] = df[gene_col]+site_id_sep+df[psite_col]
    return df

def preprocessGCT(gct_file,
                  samples_file,
                  gene_col='geneSymbol',
                  psite_col='variableSites',
                  writeToFile=False):
    df = parseGCT(gct_file, gene_col, psite_col)
    df = removeDups(df, gene_col, psite_col, samples)
    if writeToFile:
        df.to_csv('preprocessed.'+gct_file, sep='\t', index=False)
    return df
