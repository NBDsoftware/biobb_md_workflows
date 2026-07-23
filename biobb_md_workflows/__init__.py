"""biobb_md_workflows — BioBB Workflows for MD simulations with GROMACS."""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("biobb_md_workflows")
except PackageNotFoundError:  # not installed (e.g. running from a source tree)
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
