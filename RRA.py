import pandas as pd
import numpy as np
from scipy import stats


# rank and normalize guide scores
def normRank(score_vector):
  """
  ranks series, with higher values ranked higher (e.g. highest values is rank 1) and normalizes by number of values
  Returns normalized rank vector, indexed same as input vector.
  """

  normRankOrder = score_vector.rank(ascending=False) / len(score_vector.dropna())

  return normRankOrder



# for a vector of normalized ranks, return the betaScore vector, and minimum score (rho)
def betaScore_rho_p(rank_vector):
    """
    Compares each element in the vector to its corresponding value in the null distribution vector,
    using the probability mass function of the binomial distribution.
    Assigns a p-value to each element in the vector, creating the betaScore vector.
    Uses minimum betaScore as rho
    """

    rank_vector = rank_vector.dropna()
    n = len(rank_vector)

    betaScores = rank_vector.copy(deep=True)
    betaScores[0:n] = np.nan

    sorted_ranks = rank_vector.dropna().sort_values().index

    for i, k in enumerate(sorted_ranks):

        x = rank_vector[k]

        betaScore = sum([stats.binom.pmf(l, n, x, loc=0) for l in range(i, n+1)])
        betaScores[k] = betaScore

    rho = min(betaScores)

    p = min([rho*n, 1])

    return(betaScores, rho, p)
