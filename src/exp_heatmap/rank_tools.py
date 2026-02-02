"""
Rank-based transformation tools for genomic selection statistics.

This module provides functions to compute empirical rank scores from test statistics.
The rank transformation converts raw test values into genome-wide percentile ranks,
which are then transformed to -log10 scale for visualization.

Statistical Interpretation:
- Null hypothesis: Under no selection, test statistics are uniformly distributed
- Reference distribution: Genome-wide empirical distribution of the test statistic
- Transformation: rank_score = -log10(rank / total_variants)
- Tie handling: Equal values receive the same rank (standard competition ranking)

Note: These are empirical rank scores, not classical p-values. They represent the
relative extremity of each variant's test statistic compared to the genome-wide
distribution, suitable for identifying outliers but not for formal hypothesis testing.
"""

import numpy as np

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
        r = r + 1  # increase rank by one (create new rank)
        if values[i] != values[i - 1]:
            ranks.append(r)  # give new rank
        else:
            ranks.append(ranks[i - 1])  # give the same rank as before

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


# Backward compatibility aliases
compute_rank_p_vals = compute_empirical_rank_scores
compute_log_10_p_vals = compute_log10_rank_scores


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

    # get rid of results with nan (there should't be too many of them)
    test_data.dropna(axis=0, inplace=True)

    # Sort the data, so that the ranks can be given
    test_data.sort_values(by=test_name, inplace=True, ascending=top_lowest)

    test_results = test_data[test_name].values

    ranks = compute_ranks(test_results)
    test_data["rank"] = ranks

    rank_scores = compute_empirical_rank_scores(ranks)
    test_data["empirical_rank_score"] = rank_scores

    log10_rank_scores = compute_log10_rank_scores(rank_scores)
    test_data["-log10_rank_score"] = log10_rank_scores

    return test_data
