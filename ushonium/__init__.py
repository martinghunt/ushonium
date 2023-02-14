from pkg_resources import get_distribution

try:
    __version__ = get_distribution("ushonium").version
except:
    __version__ = "local"


__all__ = [
    "mafft",
    "pipeline",
    "utils",
]

from ushonium import *
