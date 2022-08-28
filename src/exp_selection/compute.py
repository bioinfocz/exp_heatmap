import exp_selection.tests as tests
import exp_selection.utils as utils


def compute(zarr_dir: str, panel_file: str, output_dir: str):
    utils.check_path_or_exit(zarr_dir)
    utils.check_path_or_exit(panel_file)

    tests.run(zarr_dir=zarr_dir, panel_file=panel_file, output_dir=output_dir)

    print()
    print(f"ZARR dir: {zarr_dir}")
    print(f"Panel file: {panel_file}")
    print(f"Output dir: {output_dir}")
