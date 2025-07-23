from exp_heatmap import tests, utils
import sys

def compute(zarr_dir: str, panel_file: str, output_dir: str, test: str, chunked: bool = False):
    utils.check_path_or_exit(zarr_dir)
    utils.check_path_or_exit(panel_file)
    try:
        tests.run(zarr_dir=zarr_dir, panel_file=panel_file, output_dir=output_dir, test=test, chunked=chunked)
    except KeyboardInterrupt:
        print("")
        sys.exit(1)

    print()
    print(f"ZARR dir: {zarr_dir}")
    print(f"Panel file: {panel_file}")
    print(f"Output dir: {output_dir}")
