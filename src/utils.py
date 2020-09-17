import re
import os
import sys


def check_path_or_exit(path):
    if not os.path.exists(path):
        print("Path {} does not exist or it is unaccessible".format(path))
        sys.exit(1)


def name_from_path(path):
    path = re.sub(r"\/+$", "", path)
    path = os.path.basename(path)
    path = re.sub(r"\.vcf$", "", path)
    path = re.sub(r"\.gz$", "", path)
    path = re.sub(r"\.zarr$", "", path)

    return path


def name_with_path(path):
    return os.path.join(os.path.dirname(path), name_from_path(path))
