import sys
import allel

import exp_heatmap.utils as utils


def prepare(recode_file: str, zarr_dir: str):
    utils.check_path_or_exit(recode_file)

    try:
        allel.vcf_to_zarr(recode_file, zarr_dir, fields="*", log=sys.stdout)
    except Exception as e:
        print(f"Error converting VCF to ZARR: {e}")
        sys.exit(1)
    
    print()
    print(f"Recoded VCF: {recode_file}")
    print(f"ZARR dir: {zarr_dir}")
