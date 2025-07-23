import sys
import allel
import zarr

from exp_heatmap import utils


def prepare(recode_file: str, zarr_dir: str) -> None:
    """
    Convert VCF file to ZARR array.
    
    Requirements:
        - zarr version < 3.0.0
    
    Args:
        - recode_file: Path to the input VCF file (recoded, SNPs only)
        - zarr_dir: Path where the ZARR directory will be created
        
    Raises:
        - SystemExit: If zarr version >= 3.0.0 is detected, if input file doesn't exist, or if conversion fails
    """
    # Check zarr version compatibility
    zarr_version = zarr.__version__
    if int(zarr_version.split('.')[0]) >= 3:
        print(f"Error: Unsupported zarr version: {zarr_version}")
        print("Please downgrade to zarr version < 3.0.0:")
        print()
        print("  pip install 'zarr<3.0.0'")
        print()
        sys.exit(1)
    
    # Check if input file exists
    utils.check_path_or_exit(recode_file)

    # Convert VCF file to ZARR array
    try:
        allel.vcf_to_zarr(recode_file, zarr_dir, fields="*", log=sys.stdout)
    except KeyboardInterrupt:
        print("")
        sys.exit(1)
    except Exception as e:
        print(f"Error converting VCF to ZARR: {e}")
        sys.exit(1)
    
    print()
    print(f"Recoded VCF: {recode_file}")
    print(f"ZARR dir: {zarr_dir}")
