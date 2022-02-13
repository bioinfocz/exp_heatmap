import exp_selection.xpehh as xpehh
import exp_selection.utils as utils


def compute(zarr_dir: str, panel_file: str, xpehh_dir: str):
    utils.check_path_or_exit(zarr_dir)
    utils.check_path_or_exit(panel_file)

    xpehh.run(zarr_dir, panel_file, xpehh_dir)

    print()
    print(f"ZARR dir: {zarr_dir}")
    print(f"Panel file: {panel_file}")
    print(f"Xpehh dir: {xpehh_dir}")
