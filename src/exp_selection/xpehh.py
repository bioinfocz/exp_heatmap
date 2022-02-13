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


def run(zarr_dir: str, panel_file: str, xpehh_dir: str):
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
    df = pd.DataFrame({"variant_pos": positions})
    df.insert(0, "name", name)

    results = []  # it will hold xpehh results of al pop pairing for given chromosome
    masks = []

    for pair in pop_pairs:
        ht_pop1 = xp_utils.get_haplotypes(gt, panel, pair[0])
        ht_pop2 = xp_utils.get_haplotypes(gt, panel, pair[1])

        print("computing XPEHH for pair " + pair[0] + " " + pair[1])
        print(
            "dimensions of haplotype data for pop "
            + pair[0]
            + ": "
            + " ".join(map(str, ht_pop1.shape))
        )
        print(
            "dimensions of haplotype data for pop "
            + pair[1]
            + ": "
            + " ".join(map(str, ht_pop2.shape))
        )
        print("dimensions of positions: " + str(len(positions)))

        result = allel.xpehh(
            h1=ht_pop1,
            h2=ht_pop2,
            pos=positions,
            map_pos=None,
            min_ehh=0.05,
            include_edges=False,
            gap_scale=20000,
            max_gap=200000,
            is_accessible=None,
            use_threads=True,
        )
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

    os.makedirs(xpehh_dir, exist_ok=True)

    for pair, res in zip(pop_pairs, results):
        result_path = os.path.join(xpehh_dir, "_".join(pair) + ".tsv")

        # add results to the dataframe with coordinates
        df["xpehh"] = res

        # Compute ascending and descending log10 rank p-values
        for order_bool in [True, False]:
            df.sort_values(by="xpehh", inplace=True, ascending=order_bool)
            test_results = df["xpehh"].values
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
