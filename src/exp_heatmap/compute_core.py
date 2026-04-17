import os

import allel
import numpy as np
import pandas as pd
import zarr

from exp_heatmap import rank_tools, utils, xp_utils
from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


def _build_output_dataframe(base_df: pd.DataFrame, test: str, result_data) -> pd.DataFrame:
    """
    Build one per-population-pair output table while preserving pair-specific missingness.

    Missing values are kept as NaN in the saved table so downstream plotting can reflect
    unavailable loci for a single pair instead of silently dropping the locus for all pairs.
    """
    pair_df = base_df.copy()
    pair_df[test] = result_data
    pair_df["-log10_p_value_ascending"] = rank_tools.compute_rank_scores_for_series(
        pair_df[test], ascending=True
    )
    pair_df["-log10_p_value_descending"] = rank_tools.compute_rank_scores_for_series(
        pair_df[test], ascending=False
    )
    return pair_df


def run(
    zarr_dir: str,
    panel_file: str,
    output_dir: str,
    test="xpehh",
    d_tajima_d_size: int = 13,
    chunked: bool = False,
):
    """
    Compute selection statistics for all population pairs and write one TSV per pair.

    Results now preserve pair-specific NaN values instead of applying a combined global
    mask across all population pairs.
    """
    logger.debug(f"Loading panel file: {panel_file}")
    panel = pd.read_csv(panel_file, sep="\t", usecols=["sample", "pop", "super_pop"])
    pop_pairs = xp_utils.create_pop_pairs(panel)

    logger.debug(f"Loading zarr data: {zarr_dir}")
    callset = zarr.open_group(zarr_dir, mode="r")

    gt, positions = xp_utils.filter_by_AF(callset, 0.05, chunked)
    samples = callset["samples"][:]

    xp_utils.check_sample_order(samples, panel["sample"])

    name = utils.name_from_path(zarr_dir)
    if test == "delta_tajima_d":
        base_df = pd.DataFrame({"variant_pos": positions[0::d_tajima_d_size][:-1]})
    else:
        base_df = pd.DataFrame({"variant_pos": positions})
    base_df.insert(0, "name", name)

    os.makedirs(output_dir, exist_ok=True)
    nonempty_pair_count = 0

    for pair in pop_pairs:
        if test in ["xpehh", "xpnsl"]:
            array_pop1 = xp_utils.get_haplotypes(gt, panel, pair[0])
            array_pop2 = xp_utils.get_haplotypes(gt, panel, pair[1])
        elif test in ["delta_tajima_d", "hudson_fst"]:
            array_pop1 = xp_utils.get_pop_allele_counts(gt, panel, pair[0])
            array_pop2 = xp_utils.get_pop_allele_counts(gt, panel, pair[1])
        else:
            raise ValueError(f"Unsupported test '{test}'")

        logger.debug(f"Computing {test.upper()} for pair {pair[0]} vs {pair[1]}")
        logger.debug(f"Population {pair[0]} dimensions: {' '.join(map(str, array_pop1.shape))}")
        logger.debug(f"Population {pair[1]} dimensions: {' '.join(map(str, array_pop2.shape))}")
        logger.debug(f"Number of positions: {len(positions)}")

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

        pair_df = _build_output_dataframe(base_df, test, result)
        nan_count = int(pair_df[test].isna().sum())
        valid_count = int(pair_df[test].notna().sum())
        if valid_count > 0:
            nonempty_pair_count += 1
        else:
            logger.warning(f"All computed values are NaN for pair {pair[0]} vs {pair[1]}")

        logger.debug(
            f"Pair {pair[0]} vs {pair[1]}: {valid_count} valid loci, {nan_count} NaN loci retained in output"
        )

        result_path = os.path.join(output_dir, "_".join(pair) + ".tsv")
        pair_df.to_csv(result_path, index=False, sep="\t")

    if nonempty_pair_count == 0:
        logger.warning("All population-pair outputs contain only NaN values.")
        if test == "delta_tajima_d":
            logger.warning("For delta Tajima's D, try increasing the 'd_tajima_d_size' parameter.")
