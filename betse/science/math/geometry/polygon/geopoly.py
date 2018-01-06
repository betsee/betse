#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional **convex polygon** (i.e., non-self-intersecting
polygon such that all line segments between any two points on the polygon
boundary remain strictly inside the polygon) functionality.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseMathPolygonException
from betse.lib.numpy import nparray
from betse.science.math.geometry import geopoint
from betse.util.type.types import type_check, SequenceTypes

# ....................{ EXCEPTIONS                         }....................
def die_unless_polygon(*polygons: SequenceTypes) -> None:
    '''
    Raise an exception unless all passed sequences are two-dimensional
    polygons.

    Parameters
    ----------
    polygons : tuple[SequenceTypes]
        Tuple of all sequences to be validated.

    Raises
    ----------
    BetseMathPointException
        If any such sequence is *not* a two-dimensional polygon.

    See Also
    ----------
    :func:`is_polygon`
        Further details.
    '''

    if not is_polygon(*polygons):
        for polygon in polygons:
            if not is_polygon(polygon):
                raise BetseMathPolygonException(
                    'Sequence not a two-dimensional polygon '
                    '(i.e., length < 3 or item length != 2): {!r}'.format(
                        polygon))

# ....................{ TESTERS                            }....................
@type_check
def is_polygon(*polygons: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences are **two-dimensional polygons** (i.e.,
    contain at least three two-dimensional points).

    Parameters
    ----------
    polygons : tuple[SequenceTypes]
        Tuple of all sequences to be tested.

    Returns
    ----------
    bool
        ``True`` only if these sequences are all two-dimensional polygons.

    See Also
    ----------
    :func:`geopoint.is_point`
        Further details on two-dimensional points.
    '''

    return all(
        len(polygon) >= 3 and geopoint.is_point(*polygon)
        for polygon in polygons
    )

# ....................{ ORIENTERS                          }....................
#FIXME: Unit test us up.
@type_check
def orient_counterclockwise(polygon: SequenceTypes) -> SequenceTypes:
    '''
    Positively reorient the passed polygon, returning a copy of this polygon
    whose vertices are **positively oriented** (i.e., sorted in
    counter-clockwise order).

    Parameters
    ----------
    polygon : SequenceTypes
        Two-dimensional sequence of all points defining the possibly non-convex
        two-dimensional polygon to be positively oriented such that:
        * The first dimension indexes each such point (in arbitrary order).
        * The second dimension indexes each coordinate of this point such that:
          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.
        Note that this function expects the type of this sequence to define an
        ``__init__()`` method accepting a passed iterable as its first and only
        positional argument. Unsurprisingly, all builtin sequences (e.g.,
        :class:`tuple`, :class:`list`) *and* :mod:`numpy` sequences (e.g.,
        :class:`numpy.ndarray`) satisfy this requirement.

    Returns
    ----------
    SequenceTypes
        Copy of this polygon positively oriented. The type of this sequence is
        the same as that of the original polygon.
    '''

    # If this sequence is *NOT* a polygon, raise an exception.
    die_unless_polygon(polygon)

    # Numpy array corresponding to this sequence. While polygon reorientation is
    # feasible in pure-Python, the Numpy-based approach is significantly more
    # efficient as the number of polygon edges increases.
    poly_verts = nparray.from_iterable(polygon)

    # Centre point of this polygon,
    poly_centre = poly_verts.mean(axis=0)

    # One-dimensional Numpy array indexing each vertex of this polygon such that
    # each element is the angle in radians between the positive X-axis and that
    # vertex, derived according to the classic mnemonic SOHCAHTOA:
    #
    #              opposite             (opposite)
    #              --------             (--------)
    #     tan(θ) = adjacent  -->  arctan(adjacent) = θ
    #
    # To sort vertices in a counter-clockwise rotation around the centre point
    # of this polygon rather than around the origin (i.e., the point (0, 0)),
    # each vertex is translated from the origin to this centre point *BEFORE*
    # obtaining this vertex's angle. Assuming standard orientation for a right
    # triangle sitting at the centre point of this polygon:
    #
    # * The opposite edge is the Y coordinate of the current vertex translated
    #   from the origin onto the centre point.
    # * The adjacent edge is the X coordinate of the current vertex translated
    #   from the origin onto the centre point.
    #
    # For safety, the np.arctan2() function intended exactly for this purpose
    # rather than the np.arctan() function intended for general-purpose
    # calculation is called. For further details, see:
    #     https://en.wikipedia.org/wiki/Atan2
    poly_angles = np.arctan2(
        poly_verts[:,1] - poly_centre[1],
        poly_verts[:,0] - poly_centre[0])

    # One-dimensional Numpy array indexing the index of each vertex of this
    # polygon, sorted in counter-clockwise order.
    poly_verts_sorted_index = poly_angles.argsort()

    # Two-dimensional Numpy array of the polygon to be returned, containing all
    # vertices sorted in this order.
    poly_verts_sorted = poly_verts[poly_verts_sorted_index]

    # Return a sequence of the same type as the passed polygon.
    return nparray.to_iterable(array=poly_verts_sorted, cls=type(polygon))
