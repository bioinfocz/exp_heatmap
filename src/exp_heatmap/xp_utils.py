"""
Specific utilities for genetic data analysis and population genetics calculations.
"""

import numpy as np
import allel
import sys
import pandas as pd
from itertools import combinations
from typing import List, Tuple, Any

def create_pop_pairs(panel: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Create all unique population pairs from a population panel.
    
    Args:
        panel: DataFrame with population information, must contain 'pop' column
        
    Returns:
        List of tuples containing population pair combinations
        
    Example:
        >>> panel = pd.DataFrame({'pop': ['AFR', 'EUR', 'ASN']})
        >>> create_pop_pairs(panel)
        [('AFR', 'EUR'), ('AFR', 'ASN'), ('EUR', 'ASN')]
    """
    populations = np.unique(panel["pop"].values)
    return list(combinations(populations, 2))

def get_haplotypes(gt_array: allel.GenotypeArray, panel: pd.DataFrame, pop: str) -> allel.HaplotypeArray:
    """
    Extract haplotype data for a specific population.
    
    Args:
        gt_array: Genotype array containing all samples
        panel: Population panel DataFrame with sample metadata
        pop: Population identifier to extract
        
    Returns:
        Haplotype array for the specified population
    """
    # Get indices of samples belonging to the specified population
    indices_pop = panel.index[panel["pop"] == pop]

    # Extract genotype data for the population
    gt_pop = gt_array.take(indices_pop, axis=1)

    return gt_pop.to_haplotypes()



def get_pop_allele_counts(gt: allel.GenotypeArray, panel: pd.DataFrame, pop: str) -> allel.AlleleCountsArray:
    """
    Calculate allele counts for a specific population.
    
    Args:
        gt: Genotype array containing all samples
        panel: Population panel DataFrame with sample metadata
        pop: Population identifier to calculate counts for
        
    Returns:
        Allele counts array for the specified population
    """
    # Get indices of samples belonging to the specified population
    indices_pop = panel.index[panel["pop"] == pop]

    # Extract genotype data for the population
    gt_pop = gt.take(indices_pop, axis=1)

    # Calculate allele counts for the population (input for PBS analysis)
    ac = gt_pop.count_alleles()
    return ac

def filter_by_AF(callset: Any, af_threshold: float, chunked: bool = False) -> Tuple[allel.GenotypeChunkedArray, allel.SortedIndex]:
    """
    Filter variants by alternate allele frequency threshold.
    
    Args:
        callset: Zarr callset containing variant data
        af_threshold: Minimum alternate allele frequency threshold
        
    Returns:
        Tuple of (filtered_genotypes, filtered_positions) where variants
        have alternate allele frequency > af_threshold
    """
    af = callset["variants/AF"][:]
    loc_variant_selection = af[:, 0] > af_threshold
    gt_zarr = callset["calldata/GT"]

    # Load the genotype as chunked array for memory efficiency
    if chunked:
        print("Using chunked array to avoid memory exhaustion")
        gt = allel.GenotypeChunkedArray(gt_zarr)
    else:
        print("Attempting to load full array into memory")
        gt = allel.GenotypeArray(gt_zarr)

    # Filter variants above the threshold
    gt_filtered = gt.compress(loc_variant_selection, axis=0)
    
    # Filter positions correspondingly
    positions = allel.SortedIndex(callset["variants/POS"])
    positions_filtered = positions.compress(loc_variant_selection, axis=0)

    return gt_filtered, positions_filtered

def check_sample_order(zarr_samples, panel_samples) -> None:
    """
    Checks the sample order between zarr and panel data.
    Exits the program if samples don't match or are in different order.
    
    Args:
        zarr_samples: Array or list of sample names from zarr file
        panel_samples: Array or pandas Series of sample names from panel file
    """
    # Convert to lists for easier manipulation
    zarr_samples_list = zarr_samples.tolist() if hasattr(zarr_samples, 'tolist') else list(zarr_samples)
    panel_samples_list = panel_samples.tolist() if hasattr(panel_samples, 'tolist') else list(panel_samples)
    
    print(f"Zarr samples: {len(zarr_samples_list)}")
    print(f"Panel samples: {len(panel_samples_list)}")
    
    # Verify that both datasets have the same number of samples
    if len(zarr_samples_list) != len(panel_samples_list):
        print(f"\nERROR: Different number of samples!")
        print(f"Zarr has {len(zarr_samples_list)} samples, panel has {len(panel_samples_list)} samples")
        sys.exit(1)
    
    # Check if samples are the same (regardless of order)
    zarr_set = set(zarr_samples_list)
    panel_set = set(panel_samples_list)
    
    if zarr_set != panel_set:
        missing_in_panel = zarr_set - panel_set
        missing_in_zarr = panel_set - zarr_set
        
        print(f"\nERROR: Different samples!")
        if missing_in_panel:
            print(f"Samples in zarr but not in panel: {missing_in_panel}")
        if missing_in_zarr:
            print(f"Samples in panel but not in zarr: {missing_in_zarr}")
        sys.exit(1)
    
    # Verify that sample order matches between datasets
    order_matches = True
    mismatches = []
    
    for i, (zarr_sample, panel_sample) in enumerate(zip(zarr_samples_list, panel_samples_list)):
        if zarr_sample != panel_sample:
            order_matches = False
            mismatches.append((i, zarr_sample, panel_sample))
    
    if order_matches:
        print(f"\nSUCCESS: Sample order is identical!")
        print(f"All {len(zarr_samples_list)} samples are in the same order.")
    else:
        print(f"\nWARNING: Sample order differs!")
        print(f"Found {len(mismatches)} mismatches:")
        
        # Display the first 10 mismatched positions for debugging
        for i, zarr_sample, panel_sample in mismatches[:10]:
            print(f"  Position {i+1}: Zarr='{zarr_sample}' vs Panel='{panel_sample}'")
        
        if len(mismatches) > 10:
            print(f"  ... and {len(mismatches) - 10} more mismatches")
        
        print(f"\nFirst 10 samples comparison:")
        print("Position\tZarr\t\tPanel")
        for i in range(min(10, len(zarr_samples_list))):
            match_status = "✓" if zarr_samples_list[i] == panel_samples_list[i] else "✗"
            print(f"{i+1}\t{zarr_samples_list[i]}\t{panel_samples_list[i]}\t{match_status}")
        
        print("\nOrder of samples in panel file does not match order of samples in zarr.")
        print("It is possible that you are using wrong panel file path e.g. from different phase than your variant data.")
        sys.exit(1)