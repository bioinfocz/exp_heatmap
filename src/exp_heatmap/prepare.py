import sys
import allel
import zarr

import exp_heatmap.utils as utils


def prepare(recode_file: str, zarr_dir: str):
    zarr_version = zarr.__version__
    if zarr_version.startswith('3.') or int(zarr_version.split('.')[0]) >= 3:
        print(f"Error: zarr version {zarr_version} is not supported.")
        print("Please downgrade to zarr version < 3.0.0:")
        print("  pip install 'zarr<3.0.0'")
        sys.exit(1)
    
    utils.check_path_or_exit(recode_file)

    try:
        allel.vcf_to_zarr(recode_file, zarr_dir, fields="*", log=sys.stdout)
    except Exception as e:
        print(f"Error converting VCF to ZARR: {e}")
        sys.exit(1)
    
    print()
    print(f"Recoded VCF: {recode_file}")
    print(f"ZARR dir: {zarr_dir}")
