__all__ = [
    "connections",
    "cvterms",
    "databases",
    "db_access",
    "entities",
    "io",
    "set_up",
    "util",
]

from importlib import import_module


def __getattr__(name):
    if name in __all__:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted([*globals(), *__all__])
