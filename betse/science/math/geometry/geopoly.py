#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional polygon functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, SequenceTypes

# ....................{ TESTERS                            }....................
#FIXME: Implement us up. To do so:
#* Choose any three adjacent polygon vertices -- say, "polygon[0]",
#  "polygon[1]", and "polygon[3]".
#* For brevity, call these points A, B, and C.
#* Then return True if and only if:
#  (Bx - Ax)(Cy - Ay) - (Cx - Ax)(By - Ay) > 0, which is equivalent to:
#  (Bx - Ax)(Cy - Ay) > (Cx - Ax)(By - Ay)
#For a partial derivation, see:
#    https://en.wikipedia.org/wiki/Curve_orientation#Practical_considerations
#FIXME: Unit test us up.
@type_check
def is_convex_orientation_counterclockwise(polygon: SequenceTypes) -> bool:
    '''
    ``True`` only if the passed convex polygon is **oriented positively** (i.e.,
    the vertices of this polygon are specified in counter-clockwise order).

    Parameters
    ----------
    polygon : SequenceTypes
        Two-dimensional sequence of all points defining the convex
        two-dimensional polygon to be tested such that:
        * The first dimension indexes each such point (in arbitrary order).
        * The second dimension indexes each coordinate of this point such that:
          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.

    Returns
    ----------
    bool
        ``True`` only if this convex polygon is oriented positively.
    '''

    pass

# ....................{ CLIPPERS                           }....................
#FIXME: Unit test us up.
@type_check
def clip_convex(
    subject_polygon: SequenceTypes, clip_polygon: SequenceTypes,
) -> SequenceTypes:
    '''
    Clip the passed subject polygon to the interior of the passed clip polygon
    by retracting all points of the former residing in the exterior of the
    latter onto the boundary of the latter.

    This function implements the Sutherland-Hodgman polygon clipping algorithm.
    Unlike most brute-force implementations of this algorithm, this
    implementation does *not* require that the vertices of this clip polygon be
    presorted in either clockwise or counter-clockwise order.

    Parameters
    ----------
    subject_polygon : SequenceTypes
        Two-dimensional sequence of all points defining the **subject polygon**
        (i.e., possibly non-convex two-dimensional polygon to be clipped) such
        that:
        * The first dimension indexes each such point (in arbitrary order).
        * The second dimension indexes each coordinate of this point such that:
          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.
    clip_polygon : SequenceTypes
        Two-dimensional sequence of all points defining the **clip polygon**
        (i.e., strictly convex two-dimensional polygon to clip the subject
        polygon against) such that:
        * The first dimension indexes each such point (in strict clockwise
          order).
        * The second dimension indexes each coordinate of this point such that:
          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.

    Returns
    ----------
    SequenceTypes
        Copy of this subject polygon clipped to this clip polygon.

    See Also
    ----------
    https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
        Pseudocode outlining this algorithm.
    '''

    # Sort the vertices of this clip polygon counter-clockwise.
    clip_polygon = orient_counterclockwise(clip_polygon)

    # Clip this subject polygon to this clip polygon oriented counter-clockwise.
    return clip_convex_counterclockwise(subject_polygon, clip_polygon)


#FIXME: Implement us up.
@type_check
def clip_convex_counterclockwise(
    subject_polygon: SequenceTypes, clip_polygon: SequenceTypes,
) -> SequenceTypes:
    '''
    Clip the passed subject polygon to the interior of the passed clip polygon
    by retracting all points of the former residing in the exterior of the
    latter onto the boundary of the latter.

    This function implements the Sutherland-Hodgman polygon clipping algorithm.
    Like most brute-force implementations of this algorithm, this implementation
    requires that the vertices of this clip polygon be presorted in
    counter-clockwise order.

    See Also
    ----------
    :func:`clip`
        Pseudocode outlining this algorithm.
    '''

    #FIXME: If this clip polygon is *NOT* convex, raise an exception.

    #FIXME: Actually raise an exception here.
    if not is_convex_orientation_counterclockwise(clip_polygon):
        pass

# ....................{ ORIENTERS                          }....................
#FIXME: Implement us up. For a useful Lua algorithm to use as a starting point,
#see the following StackOverflow question:
#    https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
#FIXME: Unit test us up.
@type_check
def orient_counterclockwise(polygon: SequenceTypes) -> SequenceTypes:
    '''
    Positively orient the passed polygon.

    This function returns a copy of this polygon whose points are guaranteed to
    be **positively oriented** (i.e., sorted in counter-clockwise order).

    Parameters
    ----------
    polygon : SequenceTypes
        Two-dimensional sequence of all points defining the possibly non-convex
        two-dimensional polygon to be positively oriented such that:
        * The first dimension indexes each such point (in arbitrary order).
        * The second dimension indexes each coordinate of this point such that:
          * The first item is the X coordinate of this point.
          * The second item is the Y coordinate of this point.

    Returns
    ----------
    SequenceTypes
        Copy of this polygon positively oriented.
    '''

    pass
