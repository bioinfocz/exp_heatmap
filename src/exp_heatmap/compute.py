from exp_heatmap import tests, utils
from exp_heatmap.logging import get_logger
import sys

logger = get_logger(__name__)

def compute(zarr_dir: str, panel_file: str, output_dir: str, test: str, chunked: bool = False):
    utils.check_path_or_exit(zarr_dir)
    utils.check_path_or_exit(panel_file)
    try:
        tests.run(zarr_dir=zarr_dir, panel_file=panel_file, output_dir=output_dir, test=test, chunked=chunked)
    except KeyboardInterrupt:
        logger.info("")
        sys.exit(1)

    logger.info("Compute completed successfully")
    logger.debug(f"ZARR dir: {zarr_dir}")
    logger.debug(f"Panel file: {panel_file}")
    logger.debug(f"Output dir: {output_dir}")
