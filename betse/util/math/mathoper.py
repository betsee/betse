#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Project-wide **mathematical operators** (i.e., low-level callables performing
core tensor operations *not* already performed out-of-the-box by either NumPy or
SciPy themselves).
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.typehints import NDArrayNdim1, NDArrayNdim2

# ....................{ INTERPOLATORS                      }....................
def cross2d(arr1: NDArrayNdim2, arr2: NDArrayNdim2) -> NDArrayNdim1:
    '''
    One-dimensional cross product of two two-dimensional input arrays.

    This function literally computes the following one-liner:

    .. code-block:: pycon

       >>> import numpy as np
       >>> arr1 = np.asarray([[1, 2], [1, 2]])
       >>> arr2 = np.asarray([[4, 5], [4, 5]])
       >>> cross2d = arr1[..., 0] * arr2[..., 1] - arr1[..., 1] * arr2[..., 0]
       >>> cross2d
       [3 -3]

    Caveats
    -------
    This function should always be called in lieu of the official
    :func:`numpy.linalg.cross` function, which no longer supports this
    two-dimensional use case for inexplicable and frankly unjustifiable reasons:

        DeprecationWarning: Arrays of 2-dimensional vectors are deprecated.
        Use arrays of 3-dimensional vectors instead. (deprecated in NumPy 2.0)

    Parameters
    ----------
    arr1: NDArrayNdim2
        First 2-dimensional input array to take the cross product of.
    arr2: NDArrayNdim2
        Second 2-dimensional input array to take the cross product of.

    Returns
    -------
    NDArrayNdim1
        One-dimensional cross product of these two-dimensional input arrays.
    '''

    return arr1[..., 0] * arr2[..., 1] - arr1[..., 1] * arr2[..., 0]
