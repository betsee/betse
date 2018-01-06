#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.science.math.geometry.geopolyconvex` submodule.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_clip() -> None:
    '''
    Unit test the :func:`betse.science.math.geometry.geopolyconvex.clip`
    function.
    '''

    # Defer heavyweight imports.
    from betse.science.math.geometry.polygon import geopolyconvex

    # Rectangle to clip all subject polygons by, oriented randomly rather than
    # counter-clockwise. Due to trivial limitations of the
    # betse.lib.numpy.nparray.to_iterable() function, this rectangle is defined
    # as a list of lists rather than tuple of tuples.
    clip_rectangle = [[2, 2], [2, -1], [-1, -1], [-1, 2],]

    # Parallelogram to be clipped, oriented randomly rather than
    # counter-clockwise.
    subject_parallelogram = [[0, 1], [2, 3], [3, 3], [1, 1]]

    # Assert a parallelogram clipped by this rectangle to produce a smaller
    # parallelogram of the same general shape. Since the
    # test_clip_counterclockwise() test exercises all other edge cases, this
    # single assertion suffices here.
    assert geopolyconvex.clip(
        subject_polygon=subject_parallelogram,
        clip_polygon=clip_rectangle,
    # Smaller parallelogram clipped by this rectangle, oriented
    # counter-clockwise. Due to issues in the underlying numpy.ndarray.tolist()
    # method internally called by this call, point coordinates are of differing
    # types. While odd, this does *NO* harm and is thus ignorable.
    ) == [[1.0, 2.0], [0, 1], [1, 1], [2.0, 2.0],]


def test_clip_counterclockwise() -> None:
    '''
    Unit test the
    :func:`betse.science.math.geometry.geopolyconvex.clip_counterclockwise`
    function.
    '''

    # Defer heavyweight imports.
    from betse.science.math.geometry.polygon import geopolyconvex

    # Rectangle to clip all subject polygons by, oriented counter-clockwise.
    clip_rectangle = ((2, 2), (-1, 2), (-1, -1), (2, -1))

    # Smaller rectangle inside this rectangle, oriented counter-clockwise.
    small_rectangle = ((1, 1), (0, 1), (0, 0), (1, 0))

    # Assert a parallelogram clipped by this rectangle to produce a smaller
    # parallelogram of the same general shape.
    assert geopolyconvex.clip_counterclockwise(
        # Parallelogram to be clipped, oriented counter-clockwise.
        subject_polygon=((3, 3), (2, 3), (0, 1), (1, 1)),
        clip_polygon=clip_rectangle,
    # Smaller parallelogram clipped by this rectangle, oriented
    # counter-clockwise.
    ) == ((2, 2), (1, 2), (0, 1), (1, 1))

    # Assert a small rectangle inside this rectangle to *NOT* clip this small
    # rectangle and hence produce the same small rectangle.
    assert geopolyconvex.clip_counterclockwise(
        subject_polygon=small_rectangle, clip_polygon=clip_rectangle,
    ) == small_rectangle

    # Assert a small rectangle outside this rectangle to fully clip this small
    # rectangle and hence produce... nothing. Absolutely nothing.
    assert geopolyconvex.clip_counterclockwise(
        # Small rectangle to be clipped, oriented counter-clockwise.
        subject_polygon=((1, 4), (0, 4), (0, 3), (1, 3)),
        clip_polygon=clip_rectangle,
    ) == ()
