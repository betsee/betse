#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level Numpy-based array testing and validation facilities.
'''

# ....................{ IMPORTS                           }....................
import numpy as np
from betse.exceptions import BetseSimUnstableNaNException
from betse.util.type.types import type_check  #, ClassType, IterableTypes
from numpy import ndarray

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_if_nan(array: ndarray) -> bool:
    '''
    Raise an exception if any element of any dimension of the passed Numpy
    array is a **NaN** (i.e., Not-a-Number).

    Parameters
    ----------
    array : ndarray
        Numpy array to be tested for NaN elements.

    Raises
    ----------

    See Also
    ----------
    :func:`is_nan`
        Further details.
    '''

    if is_nan(array):
        raise BetseSimUnstableNaNException(
            'Simulation instability detected '
            '(i.e., one or more NaNs in Numpy array). '
            'Consider a smaller time step, reduced gap junction radius, '
            'and/or reduced pump rate coefficients.'
        )

# ....................{ TESTERS                           }....................
@type_check
def is_nan(array: ndarray) -> bool:
    '''
    ``True`` only if any element of any dimension of the passed Numpy array is
    a **NaN** (i.e., Not-a-Number).

    Assuming the default Numpy configuration, Numpy produces NaN elements in
    all of the following numerical cases:

    * Division by zero (e.g., ``0/0``).
    * Division of infinity by infinity (e.g., ``inf/inf``).
    * Multplication of infinity by zero (e.g., ``0 * inf``, ``inf * 0``).
    * Subtraction of infinity by infinity (e.g., ``inf - inf``,
      ``(-inf) - (-inf)``).
    * Production of a complex number (e.g., ``np.sqrt(x)`` when ``x < 0``).
    * Performing the floating-point remainder of either:
      * A dividend that is infinite (e.g., ``np.fmod(x, y)`` when x is
        ``inf``).
      * A divisor that is zero (e.g., ``np.fmod(x, y)`` when y is 0).

    Caveats
    ----------
    **This high-level function should always be called in lieue of the
    low-level :func:`np.isnan` function,** which is space- and time-consumptive
    by compare to this function. Unlike this function, the :func:`np.isnan`
    returns a boolean array of the same shape as the passed array. In the
    common case of a two- or three-dimensional simulation matrix, passing this
    arrays to :func:`np.isnan` commonly produces a prohibitively large boolean
    array.

    This function circumvents such inefficiency by (in order):

    #. Reducing the passed array to a scalar value. If the local Numpy
       installation is linked against an optimized BLAS implementation (e.g.,
       ATLAS, OpenBLAS), this reduction is additionally parallelized across all
       available cores.
    #. Passing this scalar value to :func:`np.isnan`, producing an optimally
       space- and time-efficient scalar boolean value, which is returned as is.

    Parameters
    ----------
    array : ndarray
        Numpy array to be tested for NaN elements.

    Returns
    ----------
    bool
        ``True`` only if any element of any dimension of this array is a NaN.

    See Also
    ----------
    https://stackoverflow.com/a/6739580/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Scalar value reduced from this array. Since Numpy propagates all NaNs in
    # mathematical operations, this value is guaranteed to be either:
    #
    # * If any element of any dimension of this array is a NaN, NaN.
    # * Else, an arbitrary non-NaN value.
    #
    # To parallelize this operation across all available cores under optimized
    # BLAS implementations, this scalar is produced with the BLAS-parallelized
    # dot product operator rather than the min() or sum() functions.
    array_scalar = np.dot(array, array)

    # Return true only if this scalar value is a NaN.
    return np.isnan(array_scalar)
