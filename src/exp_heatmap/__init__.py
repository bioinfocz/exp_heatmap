# Expose main modules
from .prepare import prepare
from .compute import compute
from .plot import plot

# Expose utility modules
from . import xp_utils
from . import utils
from . import rank_tools
from . import benchmark

# Interactive module (optional, requires plotly)
try:
    from . import interactive
except ImportError:
    pass  # plotly not installed

__version__ = "1.2.0"
