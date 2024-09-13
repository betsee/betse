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
)
from beartype.vale import (
    IsAttr,
    IsEqual,
)
from numpy import ndarray

# ....................{ HINTS ~ lib : numpy                }....................
NDArrayNdim1 = Annotated[ndarray, IsAttr['ndim', IsEqual[1]]]
'''
PEP-compliant type hint matching a 1-dimensional NumPy array.
'''


NDArrayNdim2 = Annotated[ndarray, IsAttr['ndim', IsEqual[2]]]
'''
PEP-compliant type hint matching a 2-dimensional NumPy array.
'''
