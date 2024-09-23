#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Project-wide **type hints** (i.e., PEP-compliant annotations of general-purpose
interest throughout the codebase).
'''

# ....................{ IMPORTS                            }....................
from beartype.typing import (
    Annotated,
    ContextManager,
    Generator,
)
from beartype.vale import (
    IsAttr,
    IsEqual,
)
from numpy import ndarray

# ....................{ HINTS                              }....................
ContextManagerOrGenerator = ContextManager | Generator
'''
PEP-compliant type hint matching either a context manager *or* generator.
'''

# ....................{ HINTS ~ lib : numpy : 1d           }....................
NDArrayNdim1 = Annotated[ndarray, IsAttr['ndim', IsEqual[1]]]
'''
PEP-compliant type hint matching a 1-dimensional NumPy array.
'''


NDArrayNdim1Size2 = Annotated[
    ndarray,
    IsAttr['ndim', IsEqual[1]],
    IsAttr['size', IsEqual[2]],
]
'''
PEP-compliant type hint matching a 1-dimensional NumPy **2-array** (i.e., array
containing exactly two numbers, typically comprising the X and Y coordinates of
some 2-dimensional point).
'''

# ....................{ HINTS ~ lib : numpy : 2d           }....................
NDArrayNdim2 = Annotated[ndarray, IsAttr['ndim', IsEqual[2]]]
'''
PEP-compliant type hint matching a 2-dimensional NumPy array.
'''
