import allel
import zarr
import argparse
import glob
import os
from datetime import datetime
import numpy as np
import pandas as pd
import sys
import time
import json

import src.preprocessing as preprocessing
import src.xp_utils as xp_utils


def run(zarr_dir, panel_file):
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

    name = os.path.basename(zarr_dir.replace(".vcf", "").replace(".gz", ""))
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

    print("Applying NaN mask for all results")
    print("Number of results removed from each file: {}".format(num_masked))

    #  try:
   #  with open("a.json", "w") as f:
   #      for a in results:
   #          f.write(json.dumps(list(a)))

    #  except:
    #      pass

    for pair, res in zip(pop_pairs, results):
        # create direcotry for each pop pair
        pops_res_dir = "result/" + pair[0] + "_" + pair[1] + "/"
        if not os.path.exists(pops_res_dir):
            os.makedirs(pops_res_dir)

        result_path = pops_res_dir + name + ".xpehh.tsv"

        # add results to the dataframe with coordinates
        df["xpehh"] = res

        # save only the part of dataframe without nan values
        df[nan_mask].to_csv(result_path, index=False, sep="\t")

        # UPDATE LOG
        print("Resuts saved into: " + result_path)
