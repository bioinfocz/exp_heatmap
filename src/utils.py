import os
import sys


def check_path_or_exit(path):
    if not os.path.exists(path):
        print("Path {} does not exist or it is unaccessible".format(path))
        sys.exit(1)
