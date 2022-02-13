import sys
import allel

import exp_selection.utils as utils


def prepare(recode_file: str, zarr_dir: str):
    utils.check_path_or_exit(recode_file)

    allel.vcf_to_zarr(recode_file, zarr_dir, fields="*", log=sys.stdout)

    print()
    print(f"Recoded VCF: {recode_file}")
    print(f"ZARR dir: {zarr_dir}")
