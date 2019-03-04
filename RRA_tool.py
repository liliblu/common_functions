import numpy as np
from numpy import array, empty
import pandas as pd
import argparse
from scipy import stats


def fileToList(group_list):
    with open(group_list, 'r') as fh:
        return [line.strip() for line in fh.readlines()]

def cleanDF(df, sample_columns):
    '''
    Convert string nans to np.nan and string numbers to floats.
    '''
    df = df.replace(['na', 'NaN', 'Na', 'nan', 'NA', 'NAN'], np.nan)
    df[sample_columns] = df[sample_columns].astype(float)

    return df

def betaScore_rho_p(rank_vector):
    """
    Compares each element in the vector to its corresponding value in the null distribution vector,
    using the probability mass function of the binomial distribution.
    Assigns a p-value to each element in the vector, creating the betaScore vector.
    Uses minimum betaScore as rho
    """

    rank_vector = rank_vector.dropna()
    if len(rank_vector) == 0:
        return np.nan, np.nan, np.nan
    n = len(rank_vector)
    betaScores = pd.Series(empty(n))

    sorted_ranks = rank_vector.sort_values().index
    for i, k in enumerate(sorted_ranks):
        x = rank_vector[k]
        betaScore = sum([stats.binom.pmf(l, n, x, loc=0) for l in range(i, n+1)])
        betaScores[k] = betaScore

    rho = min(betaScores)
    p = min([rho*n, 1])
    return (betaScores, rho, p)

def FDR_bh(pvalues):
    pvalues = array(pvalues)
    n = sum(~np.isnan(pvalues))
    new_pvalues = empty(len(pvalues))

    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values = [x for x in values if ~np.isnan(x[0])]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in xrange(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]

    new_pvalues[np.isnan(pvalues)] = np.nan

    for i, val in enumerate(new_pvalues):
        if ~np.isnan(val):
            new_pvalues[i] = min(1, val)

    return new_pvalues


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--input_df', type=str)
    parser.add_argument('--columnsInGroup', type=str)
    parser.add_argument('--columnsOutOfGroup',type=str)
    parser.add_argument('--enrichedInHighOrLowValues', default='high', choices=['high', 'low'], type=str)
    parser.add_argument('--fdrCutOff', default=0.05, type=float)
    parser.add_argument('--outputPrefix', default='RRA', type=str)

    args = parser.parse_args()

    sample_data = pd.read_csv(args.input_df, sep='\t')
    columnsInGroup = [x for x in fileToList(args.columnsInGroup) if x in sample_data.columns]
    columnsOutOfGroup = [x for x in fileToList(args.columnsOutOfGroup) if x in sample_data.columns]
    sample_columns = columnsInGroup + columnsOutOfGroup

    enrichedInHighOrLowValues = args.enrichedInHighOrLowValues
    fdr_cut_off = args.fdrCutOff
    outputPrefix = args.outputPrefix

    sample_data = cleanDF(sample_data, sample_columns)
    ascending = enrichedInHighOrLowValues == 'low'
    sample_data[sample_columns] = sample_data[sample_columns].rank(ascending=ascending, pct=True, axis=1)

    sample_data['pval'] = sample_data[columnsInGroup].apply((lambda x: betaScore_rho_p(x)[2]), axis=1)

    sample_data['FDR'] = FDR_bh(sample_data['pval'])
    annotation_cols = [x for x in sample_data.columns if x not in sample_columns]
    sample_data[annotation_cols].to_csv('%s.RRA_results.txt' %outputPrefix, sep='\t', index=False)

    print('RRA analysis complete. Results are in %s.RRA_results.txt' %outputPrefix)
