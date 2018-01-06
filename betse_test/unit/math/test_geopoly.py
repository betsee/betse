#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.science.math.geometry.geopoly` submodule.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_orient_counterclockwise() -> None:
    '''
    Unit test the
    :func:`betse.science.math.geometry.geopoly.orient_counterclockwise`
    function.
    '''

    # Defer heavyweight imports.
    from betse.science.math.geometry.polygon import geopoly

    # Rectangle oriented randomly rather than counter-clockwise. Due to trivial
    # limitations of the betse.lib.numpy.nparray.to_iterable() function, this
    # rectangle is defined as a list of lists rather than tuple of tuples.
    rectangle_unoriented = [[2, 2], [2, -1], [-1, -1], [-1, 2],]

    # Rectangle oriented counter-clockwise.
    rectangle_oriented = geopoly.orient_counterclockwise(rectangle_unoriented)

    # Assert this rectangle to be oriented counter-clockwise.
    assert rectangle_oriented == [[-1, -1], [2, -1], [2, 2], [-1, 2],]

    # Assert that attempting to reorient this rectangle counter-clockwise is a
    # noop (i.e., produces the same rectangle).
    assert geopoly.orient_counterclockwise(rectangle_oriented) == (
        rectangle_oriented)
