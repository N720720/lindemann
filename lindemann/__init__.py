# type: ignore[attr-defined]
"""lindemann is a Python package to calculate the Lindemann index of a LAMMPS trajectory, as well as the progression of the Lindemann index per frame or per atom and frame of temperature ramps for phase transition analysis."""


try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:  # pragma: no cover
    from importlib_metadata import version, PackageNotFoundError


try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
