"""
Rank-based transformation tools for genomic selection statistics.

These helpers convert raw test values into empirical genome-wide extremeness scores
for visualization. The resulting values are useful for ranking and highlighting
outliers, but they are not calibrated inferential p-values.

Statistical interpretation:
- Reference distribution: the empirical genome-wide distribution of the chosen statistic
- Transformation: rank_score = -log10(rank / total_valid_variants)
- Tie handling: equal values receive the same rank (standard competition ranking)

These scores are suitable for visualization and prioritization, not formal hypothesis
testing under an explicit null model.
"""

import numpy as np
import pandas as pd

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


def compute_ranks(values):
    """
    Compute ranks for a sorted array of values using competition ranking.
    
    Tied values receive the same rank (the minimum rank in the tie group).
    
    Args:
        values: Sorted array of numeric values
        
    Returns:
        List of ranks corresponding to each value
    """
    i = 1
    r = 1
    ranks = [1]

    while i < len(values):
        r = r + 1
        if values[i] != values[i - 1]:
            ranks.append(r)
        else:
            ranks.append(ranks[i - 1])

        i = i + 1

    if len(ranks) == len(values):
        return ranks
    else:
        logger.error("number of ranks does not equal number of values")


def compute_empirical_rank_scores(ranks):
    """
    Convert ranks to empirical rank scores (percentile ranks).
    
    The rank score represents the fraction of variants with equal or more
    extreme values. This is NOT a classical p-value but an empirical
    genome-wide percentile.
    
    Args:
        ranks: List of integer ranks
        
    Returns:
        List of rank scores (values between 0 and 1)
    """
    ranks_size = len(ranks)
    return [rank / ranks_size for rank in ranks]


def compute_log10_rank_scores(rank_scores):
    """
    Transform rank scores to -log10 scale for visualization.
    
    Higher values indicate more extreme (potentially selected) variants.
    
    Args:
        rank_scores: List of rank scores (0-1 range)
        
    Returns:
        Numpy array of -log10 transformed scores, rounded to 3 decimals
    """
    log_scores = [np.log10(score) * -1 for score in rank_scores]
    return np.round(log_scores, 3)


def compute_rank_scores_for_series(values, ascending):
    """
    Compute -log10 empirical rank scores for one vector of test values.

    Missing values are preserved as NaN so callers can keep pair-specific missingness
    in downstream outputs and plots.
    """
    series = pd.Series(values, copy=True)
    valid_values = series.dropna()
    if valid_values.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)

    sorted_values = valid_values.sort_values(ascending=ascending)
    ranks = compute_ranks(sorted_values.values)
    rank_scores = compute_empirical_rank_scores(ranks)
    log_scores = compute_log10_rank_scores(rank_scores)

    output = pd.Series(np.nan, index=series.index, dtype=float)
    output.loc[sorted_values.index] = log_scores
    return output


def rank_across_genome(test_data, top_lowest):
    """
    Compute genome-wide empirical rank scores for test statistics.
    
    This function ranks all variants by their test statistic value and
    computes -log10 transformed rank scores for visualization.
    
    Args:
        test_data: DataFrame with test results (last column contains test values)
        top_lowest: If True, sort ascending (lowest values ranked first);
                   if False, sort descending (highest values ranked first)
                   
    Returns:
        DataFrame with added columns: 'rank', 'empirical_rank_score', '-log10_rank_score'
    """
    # The test files always have the test values saved in the last column
    test_name = test_data.columns[-1]

    valid_data = test_data.dropna(axis=0).copy()
    valid_data.sort_values(by=test_name, inplace=True, ascending=top_lowest)

    test_results = valid_data[test_name].values
    ranks = compute_ranks(test_results)
    valid_data["rank"] = ranks
    rank_scores = compute_empirical_rank_scores(ranks)
    valid_data["empirical_rank_score"] = rank_scores
    valid_data["-log10_rank_score"] = compute_log10_rank_scores(rank_scores)

    return valid_data
