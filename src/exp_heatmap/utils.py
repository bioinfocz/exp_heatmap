import re
import os
import sys


def check_path_or_exit(path: str):
    """
    Checks if the given path exists and is accessible.
    
    Prints an error message and exits the program with status code 1
    if the path does not exist or is not accessible.
    
    Args:
        path (str): The file or directory path to check
    """
    if not os.path.exists(path):
        print("Path {} does not exist or it is unaccessible".format(path))
        sys.exit(1)


def name_from_path(path: str):
    """
    Extracts a clean name from a file path by removing trailing slashes and common extensions.
    
    Args:
        path (str): The file path to process
        
    Returns:
        str: The cleaned name without directory path and extensions
    """
    path = re.sub(r"\/+$", "", path)
    path = os.path.basename(path)
    path = re.sub(r"\.vcf$", "", path)
    path = re.sub(r"\.gz$", "", path)
    path = re.sub(r"\.zarr$", "", path)

    return path

def check_sample_order(zarr_samples, panel_samples):
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
    
    # Check if same number of samples
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
    
    # Check if order is the same
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
        
        # Show first 10 mismatches
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
