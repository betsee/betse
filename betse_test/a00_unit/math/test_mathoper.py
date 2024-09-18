#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.util.math.mathoper` submodule.
'''

# ....................{ TESTS                              }....................
def test_det1d() -> None:
    '''
    Unit test the :func:`betse.util.math.mathoper.det1d` function.
    '''

    # Defer heavyweight imports.
    from betse.util.math.mathoper import det1d
    from numpy import asarray
    from numpy.testing import assert_equal

    # Two-dimensional input arrays.
    arr1 = asarray([1, 2])
    arr2 = asarray([4, 5])

    # One-dimensional cross product of these arrays.
    determinant = det1d(arr1, arr2)

    # Assert this cross product is as expected.
    assert_equal(determinant, -3, strict=True)


def test_cross2d() -> None:
    '''
    Unit test the :func:`betse.util.math.mathoper.cross2d` function.
    '''

    # Defer heavyweight imports.
    from betse.util.math.mathoper import cross2d
    from numpy import asarray
    from numpy.testing import assert_equal

    # Two-dimensional input arrays.
    arr1 = asarray([[1, 2], [6, 5]])
    arr2 = asarray([[7, 8], [3, 4]])

    # One-dimensional cross product of these arrays.
    cross_product = cross2d(arr1, arr2)

    # Assert this cross product is as expected.
    assert_equal(cross_product, asarray([-6, 9]), strict=True)
