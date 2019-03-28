#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional **convex polygon** (i.e., non-self-intersecting
polygon such that all line segments between any two points on the polygon
boundary remain strictly inside the polygon) functionality.
'''

# ....................{ IMPORTS                           }....................
import math
import numpy as np
from betse.exceptions import BetseMathPolygonException
from betse.util.type.types import type_check, SequenceTypes

# ....................{ EXCEPTIONS                        }....................
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

# ....................{ TESTERS                           }....................
@type_check
def is_polygon(*polygons: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences are **two-dimensional polygons**
    (i.e., contain at least three two-dimensional points).

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

    # Avoid circular import dependencies.
    from betse.util.math.geometry import geopoint

    # Defer to the existing is_point() tester.
    return all(
        len(polygon) >= 3 and geopoint.is_point(*polygon)
        for polygon in polygons)


@type_check
def is_convex(polygon: SequenceTypes) -> bool:
    '''
    ``True`` only if the passed polygon is **strictly convex** (i.e., is
    non-self-intersecting such that all interior angles are strictly
    between zero and a straight angle).

    Parameters
    ----------
    polygon : SequenceTypes
        Two-dimensional sequence of all points defining the polygon to be
        tested such that:

        * The first dimension indexes each such point (in arbitrary order).
        * The second dimension indexes each coordinate of this point such that:

          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.

    Returns
    ----------
    bool
        ``True`` only if this polygon is strictly convex.

    Raises
    ----------
    BetseMathPointException
        If this sequence is *not* a two-dimensional polygon.

    See Also
    ----------
    https://stackoverflow.com/a/45372025/2809027
        StackOverflow answer strongly inspiring this implementation. The credit
        entirely goes to Rory Daulton for mastering this non-trivial problem.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.numeric.floats import PI, TWO_PI

    # If this sequence is *NOT* a polygon, raise an exception.
    die_unless_polygon(polygon)

    # NOTES:  1.  Algorithm: the signed changes of the direction angles
    #             from one side to the next side must be all positive or
    #             all negative, and their sum must equal plus-or-minus
    #             one full turn (2 pi radians). Also check for too few,
    #             invalid, or repeated points.
    #         2.  No check is explicitly done for zero internal angles
    #             (180 degree direction-change angle) as this is covered
    #             in other ways, including the `n < 3` check.

    # Needed for any bad points or direction changes.
    try:
        # Get starting information.
        old_x, old_y = polygon[-2]
        new_x, new_y = polygon[-1]
        new_direction = math.atan2(new_y - old_y, new_x - old_x)
        angle_sum = 0.0

        # Check each point (the side ending there, its angle) and accumulated
        # angles.
        for ndx, newpoint in enumerate(polygon):
            # Update point coordinates and side directions and check side
            # length.
            old_x, old_y, old_direction = new_x, new_y, new_direction
            new_x, new_y = newpoint
            new_direction = math.atan2(new_y - old_y, new_x - old_x)

            # If repeated consecutive points, non-convex.
            if old_x == new_x and old_y == new_y:
                return False

            # Calculate & check the normalized direction-change angle.
            angle = new_direction - old_direction

            # Make it in half-open interval (-Pi, Pi].
            if angle <= -PI:
                angle += TWO_PI
            elif angle > PI:
                angle -= TWO_PI

            # If first time through loop, initialize orientation.
            if ndx == 0:
                if angle == 0.0:
                    return False

                orientation = 1.0 if angle > 0.0 else -1.0
            # Else if other time through loop, check orientation is stable.
            elif orientation * angle <= 0.0:  # not both pos. or both neg.
                return False

            # Accumulate the direction-change angle.
            angle_sum += angle

        # Check that the total number of full turns is plus-or-minus 1.
        return abs(round(angle_sum / TWO_PI)) == 1

    #FIXME: Clearly non-ideal. Refactor the above implementation to avoid
    #raising spurious exceptions in the first place.
    except (ArithmeticError, TypeError, ValueError):
        return False  # any exception means not a proper convex polygon


#FIXME: For orthogonality, consider refactoring to accept a single "polygon"
#parameter rather than four separate vertex parameters. Luminescent tree slime!
#FIXME: Consider raising an exception unless the four passed vertices are
#oriented counterclockwise. Fearless furballs without peers!
@type_check
def is_cyclic_quad(
    A: SequenceTypes,
    B: SequenceTypes,
    C: SequenceTypes,
    D: SequenceTypes,
) -> bool:
    '''
    ``True`` only if the quadrilateral represented by the passed four vertices
    is **cyclic** (i.e., counterclockwise-oriented).

    Parameters
    ----------
    A, B, C, D : SequenceTypes
        Four vertices of the quadrilateral to be tested, assumed to be oriented
        counterclockwise.

    Returns
    ----------
    bool
        ``True`` only if this quadrilateral is cyclic.
    '''

    # Avoid circular import dependencies.
    from betse.lib.numpy import nparray, npscalar

    # Coerce these sequences to Numpy arrays for efficiency.
    A = nparray.from_iterable(A)
    B = nparray.from_iterable(B)
    C = nparray.from_iterable(C)
    D = nparray.from_iterable(D)

    # Lengths of the four edges constructed from these vertices.
    a = np.linalg.norm(B - A)
    b = np.linalg.norm(C - B)
    c = np.linalg.norm(D - C)
    d = np.linalg.norm(A - D)

    # Calculate the length of the two diagonals.
    e = np.sqrt(((a * c + b * d) * (a * d + b * c)) / (a * b + c * d))
    f = np.sqrt(((a * c + b * d) * (a * b + c * d)) / (a * d + b * c))

    # For a cyclic quad, the product between the two diagonals equals
    # the product between the two adjacent sides.
    lhs = np.round(e*f, 15)
    rhs = np.round(a*c + b*d, 15)

    # Non-standard Numpy-specific boolean encapsulating this truth value.
    test_bool = lhs == rhs

    # Coerce this into a standard Numpy-agnostic boolean for safety.
    return npscalar.to_python(test_bool)

# ....................{ ORIENTERS                         }....................
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

    # Avoid circular import dependencies.
    from betse.lib.numpy import nparray

    # If this sequence is *NOT* a polygon, raise an exception.
    die_unless_polygon(polygon)

    # Numpy array corresponding to this sequence. While polygon reorientation
    # is feasible in pure-Python, the Numpy-based approach is significantly
    # more efficient as the number of polygon edges increases.
    poly_verts = nparray.from_iterable(polygon)

    # Centre point of this polygon,
    poly_centre = poly_verts.mean(axis=0)

    # One-dimensional Numpy array indexing each vertex of this polygon such
    # that each element is the angle in radians between the positive X-axis and
    # that vertex, derived according to the classic mnemonic SOHCAHTOA:
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
