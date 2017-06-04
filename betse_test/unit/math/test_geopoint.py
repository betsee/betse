#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.science.math.geometry.geopoint` submodule.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_is_left_of_vector() -> None:
    '''
    Unit test the :func:`betse.science.math.geometry.geopoint.is_left_of_vector`
    function.
    '''

    # Defer heavyweight imports.
    from betse.science.math.geometry import geopoint

    # Head and tail points of the vector to be tested.
    vector_head_point = [0, 2]
    vector_tail_point = [0, 0]

    # Points spatially residing to the left and right of this vector to test.
    subject_point_left  = [-1, 1]
    subject_point_right = [1, -1]

    # Ensure this function detects the correct spatiality of these points.
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
