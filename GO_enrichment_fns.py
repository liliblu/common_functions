from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import pandas as pd
import commonFns
import numpy as np


def parseGOFile(fil):
    with open(fil, 'r') as fh:
        bio_process_go = {x.strip().split('\t')[0]: set(x.strip().split('\t', 1)[1].strip().split())for x in fh.readlines()}

    M = set()
    for v in bio_process_go.values():
        M.update(set(v))
    return bio_process_go, M

def calculateGOEnrichment(genes,
                        test='hypergeom',
                        background_list=None,
                        fil='/Users/lili/dropbox_lili/common_functions/GO_terms/GO_Biological_Process_2018.txt'):

    bio_process_go, M = parseGOFile(fil)
    if background_list == None:
        M = len(M)
        background_list = M
    else:
        M = len(M.intersection(set(background_list)))

    genes = set(genes)
    N = len(genes)

    results = pd.DataFrame(columns=['GO_term', 'nGO', 'nOverlap', 'genesCommon'])
    for go, lis in bio_process_go.items():
        lis = lis.intersection(background_list)
        n = len(lis)
        if n > 0:
            k = len(lis.intersection(genes))
            line = {'GO_term':[go],'nGO':[n], 'nOverlap':[k], 'genesCommon':[','.join(lis.intersection(genes))]}
            results = results.append(pd.DataFrame.from_dict(line),sort=True, ignore_index=True)
    results['pval'] = results.apply((lambda r: hypergeom.sf(k=r['nOverlap'], M=M, n=r['nGO'], N=N, loc=1)), axis=1)
    results['FDR'] = commonFns.correct_pvalues_for_multiple_testing(results['pval'])
    results['-log10FDR'] = -np.log10(results['FDR'])
    return results
