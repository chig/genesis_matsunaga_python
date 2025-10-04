"""
Public API for genesis_interface.

Only frequently used classes are exposed at the top level.
Other modules should be imported explicitly.
"""

from typing import TYPE_CHECKING

__all__ = [
    "SMolecule",
    "STrajectories",
]

_EXPORTS = {
    "SMolecule": ("python_interface.s_molecule", "SMolecule"),
    "STrajectories": ("python_interface.s_trajectories", "STrajectories"),
}

if TYPE_CHECKING:
    from .s_molecule import SMolecule
    from .s_trajectories import STrajectories

def __getattr__(name: str):
    try:
        mod_name, attr = _EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    module = __import__(mod_name, fromlist=[attr])
    obj = getattr(module, attr)
    globals()[name] = obj  # cache for future access
    return obj

def __dir__():
    return sorted(list(globals().keys()) + __all__)

