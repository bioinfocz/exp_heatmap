import allel
import zarr
import os
import numpy as np
import pandas as pd
import sys

from exp_heatmap import xp_utils, utils, rank_tools


def run(
    zarr_dir: str,
    panel_file: str,
    output_dir: str,
    test="xpehh",
    d_tajima_d_size: int = 13,
    chunked: bool = False,
):
    """
    Computes selection tests for all population pairs.

    Args:
        zarr_dir (str): Path to VCF converted to zarr format
        panel_file (str): Path to file defining populations/super-populations  
        output_dir (str): Directory where final *.tsv files will be saved
        test (str): Type of test to compute. Options: 'xpehh', 'xpnsl', 
                   'delta_tajima_d', 'hudson_fst'. Default is 'xpehh'.
        d_tajima_d_size (int): Window size for Delta Tajima's D (number of SNPs).
                              Values between 10-20 work well. Increase if getting
                              empty results. Default is 13.
    """
    print(f"Loading panel file: {panel_file}")
    panel = pd.read_csv(panel_file, sep="\t", usecols=["sample", "pop", "super_pop"])
    pop_pairs = xp_utils.create_pop_pairs(panel)

    print(f"Loading zarr data: {zarr_dir}")
    callset = zarr.open_group(zarr_dir, mode="r")

    gt, positions = xp_utils.filter_by_AF(callset, 0.05, chunked)
    samples = callset["samples"][:]
    
    # Check the sample order
    xp_utils.check_sample_order(samples, panel["sample"])

    name = utils.name_from_path(zarr_dir)

    # For delta Tajima's D, adjust positions for non-overlapping windows
    if test == "delta_tajima_d":
        df = pd.DataFrame({"variant_pos": positions[0::d_tajima_d_size][:-1]})

    else:
        df = pd.DataFrame({"variant_pos": positions})

    df.insert(0, "name", name)

    results = []
    masks = []

    for pair in pop_pairs:
        # Prepare data structures based on test type
        if test in ["xpehh", "xpnsl"]:
            array_pop1 = xp_utils.get_haplotypes(gt, panel, pair[0])
            array_pop2 = xp_utils.get_haplotypes(gt, panel, pair[1])
        elif test in ["delta_tajima_d", "hudson_fst"]:
            array_pop1 = xp_utils.get_pop_allele_counts(gt, panel, pair[0])
            array_pop2 = xp_utils.get_pop_allele_counts(gt, panel, pair[1])

        print(f"Computing {test.upper()} for pair {pair[0]} vs {pair[1]}")
        print(f"Population {pair[0]} dimensions: {' '.join(map(str, array_pop1.shape))}")
        print(f"Population {pair[1]} dimensions: {' '.join(map(str, array_pop2.shape))}")
        print(f"Number of positions: {len(positions)}")

        # Compute the selected test
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

        # Create mask for NaN values
        mask = np.isnan(result)
        results.append(result)
        masks.append(mask)

    # Create combined NaN mask - True where any population pair has NaN
    combined_nan_mask = np.logical_or.reduce(masks)
    # Convert to boolean mask for valid positions (not NaN)
    valid_positions_mask = ~combined_nan_mask
    
    num_masked = np.sum(combined_nan_mask)
    print("Applying NaN mask for all results")
    print("Number of results removed from each file: {}".format(num_masked))

    # Check if all positions are masked
    if num_masked == len(valid_positions_mask):
        print("\n=== WARNING ===")
        print("All positions have NaN results. No output will be generated.")
        if test == "delta_tajima_d":
            print("For delta Tajima's D, try increasing the 'd_tajima_d_size' parameter.")
        return

    os.makedirs(output_dir, exist_ok=True)

    for pair, result_data in zip(pop_pairs, results):
        result_path = os.path.join(output_dir, "_".join(pair) + ".tsv")

        # add results to the dataframe with coordinates
        df[test] = result_data

        # Compute ascending and descending log10 rank p-values
        for ascending in [True, False]:
            df.sort_values(by=test, inplace=True, ascending=ascending)
            test_results = df[test].values
            ranks = rank_tools.compute_ranks(test_results)
            rank_p_vals = rank_tools.compute_rank_p_vals(ranks)
            log_10_p_vals = rank_tools.compute_log_10_p_vals(rank_p_vals)

            suffix = "ascending" if ascending else "descending"
            df[f"-log10_p_value_{suffix}"] = log_10_p_vals

        df.sort_values(by="variant_pos", inplace=True, ascending=True)

        # Save only positions without NaN values
        df[valid_positions_mask].to_csv(result_path, index=False, sep="\t")

        print(f"Results saved to: {result_path}")
