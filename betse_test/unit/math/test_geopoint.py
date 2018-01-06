#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.science.math.geometry.geopoint` submodule.
'''

# ....................{ IMPORTS                            }....................
import pytest

# ....................{ TESTS                              }....................
def test_is_left_of_vector() -> None:
    '''
    Unit test the :func:`betse.science.math.geometry.geopoint.is_left_of_vector`
    function.
    '''

    # Defer heavyweight imports.
    from betse.science.math.geometry import geopoint

    # Head and tail points of the vector to be tested.
    vector_head_point = (0, 2)
    vector_tail_point = (0, 0)

    # Points spatially residing to the left and right of this vector to test.
    subject_point_left  = (-1, 1)
    subject_point_right = (1, -1)

    # Assert the expected spatiality of these points.
    assert geopoint.is_left_of_vector(
        subject_point=subject_point_left,
        vector_head_point=vector_head_point,
        vector_tail_point=vector_tail_point,
    )
    assert not geopoint.is_left_of_vector(
        subject_point=subject_point_right,
        vector_head_point=vector_head_point,
        vector_tail_point=vector_tail_point,
    )


def test_intersect_lines() -> None:
    '''
    Unit test the :func:`betse.science.math.geometry.geopoint.intersect_lines`
    function.
    '''

    # Defer heavyweight imports.
    import numpy
    from betse.exceptions import BetseMathLineException
    from betse.science.math.geometry import geopoint
    from numpy import array

    # Assert the expected intersection of two intersecting lines given tuples.
    assert geopoint.intersect_lines(
        line1_point1=( 1,  2),
        line1_point2=(-2, -1),
        line2_point1=(-2,  1),
        line2_point2=( 1, -2),
    ) == (-1, 0)

    # Assert the same intersection of the same intersecting lines given Numpy
    # arrays.
    intersection_array = geopoint.intersect_lines(
        line1_point1=array(( 1,  2)),
        line1_point2=array((-2, -1)),
        line2_point1=array((-2,  1)),
        line2_point2=array(( 1, -2)),
    )
    assert numpy.array_equal(intersection_array, array((-1, 0)))

    # Assert two collinear lines to have no unique intersection.
    with pytest.raises(BetseMathLineException):
        geopoint.intersect_lines(
            line1_point1=( 1,  2),
            line1_point2=(-2, -1),
            line2_point1=( 0,  1),
            line2_point2=(-1,  0),
        )

    # Assert two parallel lines to have no intersection at all.
    with pytest.raises(BetseMathLineException):
        geopoint.intersect_lines(
            line1_point1=( 1,  2),
            line1_point2=(-2, -1),
            line2_point1=( 1,  0),
            line2_point2=( 0, -1),
        )
