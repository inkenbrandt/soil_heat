"""
soil_heat: A Python library for soil heat flux models.
======================================================

This package provides a collection of Python functions for calculating
soil heat flux and related soil thermal properties. The implementations
are based on various models and parameterizations published in the
peer-reviewed scientific literature.

The library is organized into modules, each corresponding to a specific
publication or a set of related equations.

Available submodules:
---------------------
- `soil_heat`: Core functions and utilities.
- `gao_et_al`: Functions from Gao et al. (2017).
- `liebethal_and_folken`: Functions from Liebethal & Foken (2006).
- `wang_and_bouzeid`: Functions from Wang & Bou-Zeid (2012).
- `wang_and_yang`: Functions from Yang & Wang (2008).

All functions are designed to work with NumPy arrays for efficient,
vectorized computations.

"""

__author__ = """Paul Inkenbrandt"""
__email__ = "paulinkenbrandt@utah.gov"
__version__ = "0.1.1"

from .soil_heat import *
from .gao_et_al import *
from .liebethal_and_folken import *
from .wang_and_bouzeid import *
from .wang_and_yang import *
