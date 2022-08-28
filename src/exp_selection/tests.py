import allel
import zarr
import os
import numpy as np
import pandas as pd
import sys

import exp_selection.preprocessing as preprocessing
import exp_selection.xp_utils as xp_utils
import exp_selection.utils as utils
import exp_selection.rank_tools as rank_tools


def run(
    zarr_dir: str,
    panel_file: str,
    output_dir: str,
    test="xpehh",
    d_tajima_d_size=13,
):
    """
    Computes the selection test for all your population pairs.

    zarr_dir -- vcf turned to zarr

    panel_file -- population/super-populations defined

    output_dir -- directory where the final *.tsv files end up

    test -- what kind of test to compute, using the scikit-allele selection or pairwise distances/F-statistics tests.
            Should be one of the following: 'xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst', 'custom'

    d_tajima_d_size -- int, Delta Tajima's D size and step value (number of SNPs), values between 10 and 20 seems to work well
                       and are sufficiently large not to be too sensitive to local variation on one side, but still be able to describe
                       subtle local changes in Tajima's D.
                       If you are getting no results (empty *.tsv), try increasing this number to be sure, that there are enough SNPs
                       between the all possible population pairs to compute Tajima's D.


    """

    panel = pd.read_csv(panel_file, sep="\t", usecols=["sample", "pop", "super_pop"])
    pop_pairs = xp_utils.create_pop_pairs(panel)

    callset = zarr.open_group(zarr_dir, mode="r")

    gt, positions = preprocessing.filter_by_AF(callset, 0.05)

    samples = callset["samples"][:]
    if np.all(samples == panel["sample"].values):
        print("Order of samples ok")
    else:
        print(
            "Order of samples in panel file does not match order of samples in given zarr. It is possible that you are using wrong panel file path e.g. from different phase than you variant data comes from different phase than your data"
        )

        sys.exit(1)

    name = utils.name_from_path(zarr_dir)

    # in delta Tajima's D we need to compute the test on bunch of SNPs, the resulting table (df) will be shorter
    # than the input data
    # here, we adjust the positions according to moving non-overlapping window of d_tajima_d_size SNPs
    if test == "delta_tajima_d":
        df = pd.DataFrame({"variant_pos": positions[0::d_tajima_d_size][:-1]})

    else:
        df = pd.DataFrame({"variant_pos": positions})

    df.insert(0, "name", name)

    results = (
        []
    )  # it will hold xpehh results of al pop pairing for given chromosomeXPNSL
    masks = []

    for pair in pop_pairs:

        # prepare datastructures for population pair GTs, depending on the test
        # haplotype array
        if test in ["xpehh", "xpnsl"]:
            array_pop1 = xp_utils.get_haplotypes(gt, panel, pair[0])
            array_pop2 = xp_utils.get_haplotypes(gt, panel, pair[1])

        # allele counts array
        elif test in ["delta_tajima_d", "hudson_fst"]:
            array_pop1 = xp_utils.get_pop_alllel_counts(gt, panel, pair[0])
            array_pop2 = xp_utils.get_pop_alllel_counts(gt, panel, pair[1])

        else:
            print("test not set")
            sys.exit()

        print("computing {} for pair ".format(test.upper()) + pair[0] + " " + pair[1])
        print(
            "dimensions of haplotype data for pop "
            + pair[0]
            + ": "
            + " ".join(map(str, array_pop1.shape))
        )
        print(
            "dimensions of haplotype data for pop "
            + pair[1]
            + ": "
            + " ".join(map(str, array_pop2.shape))
        )
        print("dimensions of positions: " + str(len(positions)))

        if test == "xpehh":

            result = allel.xpehh(
                h1=array_pop1,
                h2=array_pop2,
                pos=positions,
                map_pos=None,
                min_ehh=0.05,
                include_edges=False,
                gap_scale=20000,
                max_gap=200000,
                is_accessible=None,
                use_threads=True,
            )

        elif test == "xpnsl":

            result = allel.xpnsl(
                h1=array_pop1,
                h2=array_pop2,
                use_threads=True,
            )

        elif test == "delta_tajima_d":

            result = allel.moving_delta_tajima_d(
                ac1=array_pop1,
                ac2=array_pop2,
                size=d_tajima_d_size,
                start=0,
                stop=None,
                step=d_tajima_d_size,
            )

        elif test == "hudson_fst":

            num, den = allel.hudson_fst(
                ac1=array_pop1,
                ac2=array_pop2,
            )
            result = num / den

        else:
            print("test not set")
            sys.exit()

        mask = np.isnan(result)

        results.append(result)
        masks.append(mask)

    # create the final nan mask that will mask every position where nan occured
    # in any of the pop pairing xpehh
    # initialize the mask with one of the masks in masks
    nan_mask = masks[0]

    # then compare final nan_mask with each mask to store True whenever there is True in either mask
    for m in masks:
        nan_mask = nan_mask | m

    # finally, I will negate the whole mask bc I actually want to have
    # False in places where there is NaN
    nan_mask = [not i for i in nan_mask]

    # count the number of results that will be removed from each file after masking
    num_masked = nan_mask.count(False)

    print(nan_mask)

    print("Applying NaN mask for all results")
    print("Number of results removed from each file: {}".format(num_masked))

    # stop, if all positions are masked --> no results will be saved anyway
    if num_masked == len(nan_mask):

        print()
        print("=== !!! ===")
        print(
            "All positions are masked as having NaN results. No output. If you are computing delta Tajima's D, try increasing the 'tajima_d_size' parameter."
        )

        return

    os.makedirs(output_dir, exist_ok=True)

    for pair, res in zip(pop_pairs, results):
        result_path = os.path.join(output_dir, "_".join(pair) + ".tsv")

        # add results to the dataframe with coordinates
        df[test] = res

        # Compute ascending and descending log10 rank p-values
        for order_bool in [True, False]:
            df.sort_values(by=test, inplace=True, ascending=order_bool)
            test_results = df[test].values
            ranks = rank_tools.compute_ranks(test_results)
            rank_p_vals = rank_tools.compute_rank_p_vals(ranks)
            log_10_p_vals = rank_tools.compute_log_10_p_vals(rank_p_vals)

            if order_bool:
                df["-log10_p_value_ascending"] = log_10_p_vals
            else:
                df["-log10_p_value_descending"] = log_10_p_vals

        df.sort_values(by="variant_pos", inplace=True, ascending=True)

        # save only the part of dataframe without nan values
        df[nan_mask].to_csv(result_path, index=False, sep="\t")

        # UPDATE LOG
        print("Resuts saved into: " + result_path)
