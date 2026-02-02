import sys
import allel
import zarr

from exp_heatmap import utils
from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


class LoggerWriter:
    """Adapter to make a logging.Logger behave like a file-like object for scikit-allel."""
    def __init__(self, logger_instance, level=None):
        self.logger = logger_instance
        self.level = level or logger_instance.info
    
    def write(self, message):
        # Strip trailing newlines to avoid double-spacing in logs
        message = message.rstrip('\n')
        if message:  # Only log non-empty messages
            self.level(message)
    
    def flush(self):
        # Required for file-like interface, but logging handles this
        pass


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
        logger.error(f"Unsupported zarr version: {zarr_version}")
        logger.error("Please downgrade to zarr version < 3.0.0:")
        logger.error("  pip install 'zarr<3.0.0'")
        sys.exit(1)
    
    # Check if input file exists
    utils.check_path_or_exit(recode_file)

    # Convert VCF file to ZARR array
    try:
        # Wrap logger in file-like adapter for scikit-allel compatibility
        log_writer = LoggerWriter(logger, level=logger.debug)
        allel.vcf_to_zarr(recode_file, zarr_dir, fields="*", log=log_writer)
    except KeyboardInterrupt:
        logger.info("")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error converting VCF to ZARR: {e}")
        sys.exit(1)
    
    logger.info(f"Prepare completed successfully")
    logger.debug(f"Recoded VCF: {recode_file}")
    logger.debug(f"ZARR dir: {zarr_dir}")
