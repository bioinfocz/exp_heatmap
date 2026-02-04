"""
Logging configuration module for ExP Heatmap.

This module provides centralized logging setup that separates verbose technical
output (redirected to log files) from essential user-facing status messages
(kept in stdout). Each CLI command run generates a timestamped log file.

Usage:
    from exp_heatmap.logging import setup_logging, get_logger
    
    # In CLI entry point:
    setup_logging("compute", log_dir="logs", verbose=False)
    
    # In any module:
    logger = get_logger(__name__)
    logger.info("Status message")  # Goes to console + file
    logger.debug("Technical detail")  # Goes to file only (unless verbose=True)
"""

import logging
import os
from datetime import datetime
from typing import Optional

LOGGER_NAME = "exp_heatmap"

_logging_initialized = False
_log_file_path: Optional[str] = None


def setup_logging(
    command_name: str,
    log_dir: str = "logs",
    verbose: bool = False,
    log_to_file: bool = True
) -> logging.Logger:
    """
    Initialize logging for a CLI command run.
    
    Creates a timestamped log file and configures console output.
    Console shows INFO and above by default, DEBUG if verbose=True.
    File always captures DEBUG and above.
    
    Parameters
    ----------
    command_name : str
        Name of the CLI command (e.g., 'prepare', 'compute', 'plot').
        Used in the log filename.
    log_dir : str, optional
        Directory for log files (default: 'logs' in current working directory).
    verbose : bool, optional
        If True, show DEBUG messages in console (default: False).
    log_to_file : bool, optional
        If True, write logs to file (default: True).
        
    Returns
    -------
    logging.Logger
        Configured root logger for the package.
    """
    global _logging_initialized, _log_file_path
    
    logger = logging.getLogger(LOGGER_NAME)
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    
    # Create formatters
    console_formatter = logging.Formatter(
        fmt="%(message)s"
    )
    file_formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    
    # Handlers
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    if log_to_file:
        os.makedirs(log_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = f"exp_heatmap_{command_name}_{timestamp}.log"
        _log_file_path = os.path.join(log_dir, log_filename)
        
        file_handler = logging.FileHandler(_log_file_path, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    _logging_initialized = True
    
    # Log initial setup info
    logger.debug(f"Logging initialized for command: {command_name}")
    if log_to_file and _log_file_path:
        logger.debug(f"Log file: {_log_file_path}")
    
    return logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for a specific module.
    
    If logging hasn't been set up yet, returns a basic logger that
    outputs to console only (useful for library usage without CLI).
    
    Parameters
    ----------
    name : str
        Module name, typically __name__ (e.g., 'exp_heatmap.plot').
        
    Returns
    -------
    logging.Logger
        Logger instance for the module.
    """
    if name.startswith(LOGGER_NAME):
        logger_name = name
    else:
        module_name = name.split('.')[-1] if '.' in name else name
        logger_name = f"{LOGGER_NAME}.{module_name}"
    
    logger = logging.getLogger(logger_name)
    
    # Set handlers if logging hasn't been initialized
    if not _logging_initialized and not logger.handlers:
        parent_logger = logging.getLogger(LOGGER_NAME)
        if not parent_logger.handlers:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)
            console_handler.setFormatter(logging.Formatter("%(message)s"))
            parent_logger.addHandler(console_handler)
            parent_logger.setLevel(logging.DEBUG)
    
    return logger


def get_log_file_path() -> Optional[str]:
    """
    Get the path to the current log file.
    
    Returns
    -------
    str or None
        Path to the log file, or None if file logging is not enabled.
    """
    return _log_file_path


def set_console_level(level: int) -> None:
    """
    Change the console logging level at runtime.
    
    Parameters
    ----------
    level : int
        Logging level (e.g., logging.DEBUG, logging.INFO).
    """
    logger = logging.getLogger(LOGGER_NAME)
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
            handler.setLevel(level)

