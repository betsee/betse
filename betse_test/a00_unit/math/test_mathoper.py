#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.math.mathoper` submodule.
'''

# ....................{ TESTS                              }....................
def test_cross2d() -> None:
    '''
    Unit test the :func:`betse.util.math.mathoper.cross2d` function.
    '''

    # Defer heavyweight imports.
    from betse.util.math.mathoper import cross2d
    from numpy import asarray

    # Two-dimensional input arrays.
    arr1 = asarray([[1, 2], [1, 2]])
    arr2 = asarray([[4, 5], [4, 5]])

    # One-dimensional cross product of these arrays.
    cross_product = cross2d(arr1, arr2)

    # Assert this cross product is as expected.
    assert cross_product == asarray([3, -3])
