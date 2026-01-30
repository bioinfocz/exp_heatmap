# Expose main modules
from .prepare import prepare
from .compute import compute
from .plot import plot

# Expose utility modules
from . import xp_utils
from . import utils
from . import rank_tools

#TODO: ADD THE MODULE BELOW ONLY IF IMPLEMENTED INTO PROD

from . import benchmark
from . import interactive

__version__ = "1.2.0"
