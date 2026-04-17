from importlib import import_module

__version__ = "1.3.0"

_LAZY_EXPORTS = {
    "prepare": ("prepare", "prepare"),
    "compute": ("compute", "compute"),
    "plot": ("plot", "plot"),
    "xp_utils": ("xp_utils", None),
    "utils": ("utils", None),
    "rank_tools": ("rank_tools", None),
    "exp_logging": ("logging", None),
    "setup_logging": ("logging", "setup_logging"),
    "get_logger": ("logging", "get_logger"),
    "interactive": ("interactive", None),
}

__all__ = ["__version__", *_LAZY_EXPORTS]


def __getattr__(name):
    try:
        module_name, attr_name = _LAZY_EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'") from exc

    module = import_module(f".{module_name}", __name__)
    value = module if attr_name is None else getattr(module, attr_name)
    globals()[name] = value
    return value
