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
from betse.util.type.typehints import (
    NDArrayNdim1,
    NDArrayNdim1Size2,
    NDArrayNdim2,
)
from numbers import Number

# ....................{ OPERATIONS ~ 1d                    }....................
def det1d(arr1: NDArrayNdim1Size2, arr2: NDArrayNdim1Size2) -> Number:
    '''
    Scalar determinant of two one-dimensional input **2-arrays** (i.e., arrays
    containing exactly two numbers each).

    This function literally computes the following one-liner:

    .. code-block:: pycon

       >>> import numpy as np
       >>> arr1 = np.asarray([1, 2])
       >>> arr2 = np.asarray([4, 5])
       >>> det1d = arr1[0] * arr2[1] - arr1[1] * arr2[0]
       >>> det1d
       -3

    Caveats
    -------
    This function is often (but *not* always) called in lieu of the official
    :func:`numpy.linalg.det` function, which requires that both of the passed
    2-arrays be congealed into a single two-dimensional input array. Weird!

    Parameters
    ----------
    arr1: NDArrayNdim2
        First 1-dimensional input 2-array to take the determinant of.
    arr2: NDArrayNdim2
        Second 1-dimensional input 2-array to take the determinant of.

    Returns
    -------
    Number
        Scalar determinant of these one-dimensional input 2-arrays.
    '''

    # Return us up the one-liner for great justice!
    return arr1[0] * arr2[1] - arr1[1] * arr2[0]

# ....................{ OPERATIONS ~ 2d                    }....................
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

    # Return us up the one-liner bomb.
    return arr1[..., 0] * arr2[..., 1] - arr1[..., 1] * arr2[..., 0]
