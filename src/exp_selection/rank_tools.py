import numpy as np


def compute_ranks(values):
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
        print("number of ranks does not equal number of values")


def compute_rank_p_vals(ranks):
    ranks_size = len(ranks)
    return [rank / ranks_size for rank in ranks]


def compute_log_10_p_vals(pvals):
    log_pvals = [np.log10(pval) * -1 for pval in pvals]
    return np.round(log_pvals, 3)


def rank_across_genome(test_data, top_lowest):

    # The test files always have the test values saved in the last column
    test_name = test_data.columns[-1]

    # get rid of results with nan (there should't be too many of them)
    test_data.dropna(axis=0, inplace=True)

    # Sort the data, so that the ranks can be given
    test_data.sort_values(by=test_name, inplace=True, ascending=top_lowest)

    test_results = test_data[test_name].values

    ranks = compute_ranks(test_results)
    test_data["rank"] = ranks

    rank_p_vals = compute_rank_p_vals(ranks)
    test_data["rank_p_value"] = rank_p_vals

    log_10_p_vals = compute_log_10_p_vals(rank_p_vals)
    test_data["-log10_p_value"] = log_10_p_vals

    return test_data
