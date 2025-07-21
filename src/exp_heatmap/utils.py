"""
General utilities for file path handling and name extraction.
"""

import re
import os
import sys


def check_path_or_exit(path: str):
    """
    Checks if the given path exists and is accessible.
    
    Prints an error message and exits the program with status code 1
    if the path does not exist or is not accessible.
    
    Args:
        path (str): The file or directory path to check
    """
    if not os.path.exists(path):
        print("Path {} does not exist or it is unaccessible".format(path))
        sys.exit(1)


def name_from_path(path: str):
    """
    Extracts a clean name from a file path by removing trailing slashes and common extensions.
    
    Args:
        path (str): The file path to process
        
    Returns:
        str: The cleaned name without directory path and extensions
    """
    path = re.sub(r"\/+$", "", path)
    path = os.path.basename(path)
    path = re.sub(r"\.vcf$", "", path)
    path = re.sub(r"\.gz$", "", path)
    path = re.sub(r"\.zarr$", "", path)

    return path
