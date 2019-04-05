from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import pandas as pd
import commonFns

fil = '/Users/lili/google_drive/Ruggles_lab/common_functions/GO_terms/go_biological_process_gene_map.txt'
with open(fil, 'r') as fh:
    bio_process_go = {x.split('|')[0]: set(x.split('|')[1].split(' ')) for x in fh.readlines()}

M = set()
for v in bio_process_go.values():
    M.update(set(v))
M = len(M)

def calculateGOEnrichment(genes, test='hypergeom'):

    genes = set(genes)
    N = len(genes)
    results = pd.DataFrame(columns=['GO_term', 'nGO', 'nOverlap', 'genesGO', 'genesCommon'])
    for go, lis in bio_process_go.iteritems():
        n = len(lis)
        k = len(lis.intersection(genes))
        line = {'GO_term':[go],'nGO':[n], 'nOverlap':[k],
        'genesGO':[','.join(lis)], 'genesCommon':[','.join(lis.intersection(genes))]}
        results = results.append(pd.DataFrame.from_dict(line),sort=True, ignore_index=True)
    results['pval'] = results.apply((lambda r: hypergeom.sf(k=r['nOverlap'], M=M, n=r['nGO'], N=N, loc=1)), axis=1)
    results['FDR'] = commonFns.correct_pvalues_for_multiple_testing(results['pval'])
    return results
