#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional **convex polygon** (i.e., non-self-intersecting
polygon such that all line segments between any two points on the polygon
boundary remain strictly inside the polygon) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.geometry import geopoint
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
def is_orientation_counterclockwise(polygon: SequenceTypes) -> bool:
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
def clip(
    subject_polygon: SequenceTypes, clip_polygon: SequenceTypes,
) -> SequenceTypes:
    '''
    Clip the passed possibly non-convex subject polygon to the interior of the
    passed convex clip polygon by retracting all points of the former residing
    in the exterior of the latter onto the boundary of the latter.

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
    '''

    # Sort the vertices of this clip polygon counter-clockwise.
    clip_polygon = orient_counterclockwise(clip_polygon)

    # Clip this subject polygon to this clip polygon oriented counter-clockwise.
    return clip_counterclockwise(subject_polygon, clip_polygon)


#FIXME: Unit test us up.
#FIXME: Optimize. At the least, the geopoint.is_left_of_vector() function
#should be inlined; ideally, so should geopoint.intersect_lines()
@type_check
def clip_counterclockwise(
    subject_polygon: SequenceTypes, clip_polygon: SequenceTypes,
) -> SequenceTypes:
    '''
    Clip the passed possibly non-convex subject polygon to the interior of the
    passed convex clip polygon by retracting all points of the former residing
    in the exterior of the latter onto the boundary of the latter.

    This function implements the Sutherland-Hodgman polygon clipping algorithm.
    Like most brute-force implementations of this algorithm, this implementation
    requires that the vertices of this clip polygon be presorted in
    counter-clockwise order.

    See Also
    ----------
    https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
        Pseudocode outlining this algorithm.
    '''

    #FIXME: If this clip polygon is *NOT* convex, raise an exception.

    #FIXME: Actually raise an exception here.
    # if not is_convex_orientation_counterclockwise(clip_polygon):
    #     pass

    # Current version of the clipped subject polygon to be returned, where
    # current implies the version of this subject polygon to be used in the
    # current iteration of the inner loop below.
    subj_poly_curr = subject_polygon

    # Previous version of the clipped subject polygon to be returned, where
    # previous implies the version of this subject polygon used in the prior
    # iteration of the inner loop below.
    subj_poly_prev = None

    # Current vertex of the clip polygon to clip against. The line segment
    # between this and the next such vertex is the current edge of this clip
    # polygon to clip the current edge of this subject polygon against.
    clip_vert_curr = clip_polygon[-1]

    # Next vertex of the clip polygon to clip against.
    clip_vert_next = None

    # Current vertex of the subject polygon to be clipped. The line segment
    # between this and the next such vertex is the current edge of this subject
    # polygon to be clipped against the current edge of this clip polygon.
    subj_vert_curr = None

    # Next vertex of the subject polygon to be clipped.
    subj_vert_next = None

    # True only if the current subject vertex resides inside this clip polygon.
    is_subj_vert_curr_inside = None

    # True only if the next subject vertex resides inside this clip polygon.
    is_subj_vert_next_inside = None

    # For each next vertex of the clip polygon to clip against...
    for clip_vert_next in clip_polygon:
        # Preserve the previous version of the clipped subject polygon produced
        # by the previous iteration of the inner loop below *BEFORE* destroying
        # this sequence.
        subj_poly_prev = subj_poly_curr

        # Current version of the clipped subject polygon to be returned,
        # initialized to the empty list. For each edge of the previous version
        # of the clipped subject polygon at least partially residing inside the
        # clip polygon, the next iteration of the inner loop below appends at
        # least one vertex along this edge to this sequence.
        subj_poly_curr = []

        # Current vertex of the previous version of the clipped subject polygon
        # to be clipped, initialized to the last vertex of this polygon.
        subj_vert_curr = subj_poly_prev[-1]

        # This vertex resides inside this clip polygon if and only if this
        # vertex is spatially situated to the left of the vector signifying the
        # current edge of this clip polygon oriented in the counter-clockwise
        # orientation of this clip polygon, whose:
        #
        # * Head is the next vertex of this clip polygon.
        # * Tail is the current vertex of this clip polygon.
        is_subj_vert_curr_inside = geopoint.is_left_of_vector(
            subj_point=subj_vert_curr,
            vector_head_point=clip_vert_next,
            vector_tail_point=clip_vert_curr,)

        # For each next vertex of the previous version of the clipped subject
        # polygon...
        for subj_vert_next in subj_poly_prev:
            # This vertex resides inside this clip polygon if and only if this
            # vertex is spatially situated to the left of the same vector.
            is_subj_vert_next_inside = geopoint.is_left_of_vector(
                subj_point=subj_vert_next,
                vector_head_point=clip_vert_next,
                vector_tail_point=clip_vert_curr,)

            # If the current and next vertices of the previous version of the
            # clipped subject polygon reside on different sides of the current
            # edge of the clip polygon, the current edge of the former
            # necessarily intersects the current edge of the latter. In this
            # case, retract whichever of these two vertices resides outside of
            # the clip polygon to the boundary of the clip polygon.
            if is_subj_vert_curr_inside != is_subj_vert_next_inside:
                subj_poly_curr.append(geopoint.intersect_lines(
                    line1_point1=clip_vert_curr,
                    line1_point2=clip_vert_next,
                    line2_point1=subj_vert_curr,
                    line2_point2=subj_vert_next,
                ))

            # If the next vertex of the previous version of the clipped subject
            # polygon resides inside the current edge of the clip polygon,
            # preserve this vertex. To preserve the current orientation of the
            # subject polygon, this is performed *AFTER* possibly retracting a
            # subject vertex above.
            #
            # Note that similar logic need *NOT* be performed for the current
            # vertex of the previous version of the clipped subject polygon.
            # Why? Because if this vertex also resides inside the current edge
            # of the clip polygon, this vertex either:
            #
            # * Is the last vertex, in which case the last iteration of this
            #   loop is already guaranteed to preserve this vertex.
            # * Is any other vertex, in which case the prior iteration of this
            #   loop already preserved this vertex.
            #
            # In either case, this vertex is already preserved.
            if is_subj_vert_next_inside:
                subj_poly_curr.append(subj_vert_next)

            # Rotate the current vertex of the previous version of the clipped
            # subject polygon to the next such vertex.
            subj_vert_curr              = subj_vert_next
            is_subj_vert_curr_inside = is_subj_vert_next_inside

        # Rotate the current vertex of the clip polygon to the next such vertex.
        clip_vert_curr = clip_vert_next

    # Return the current and hence final version of the clipped subject polygon.
    return subj_poly_curr

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
