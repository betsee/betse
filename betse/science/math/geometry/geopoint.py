#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional point functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMathLineException, BetseMathPointException
from betse.util.type import iterables
from betse.util.type.types import type_check, SequenceTypes

# ....................{ EXCEPTIONS                         }....................
def die_unless_point(*points: SequenceTypes) -> None:
    '''
    Raise an exception unless all passed sequences are two-dimensional points.

    Parameters
    ----------
    points : tuple[SequenceTypes]
        Tuple of all sequences to be validated.

    Raises
    ----------
    BetseMathPointException
        If any such sequence is *not* a two-dimensional point.

    See Also
    ----------
    :func:`is_point`
        Further details.
    '''

    if not is_point(*points):
        for point in points:
            if not is_point(point):
                raise BetseMathPointException(
                    'Sequence not a two-dimensional point '
                    '(i.e., length != 2): {!r}'.format(point))

# ....................{ TESTERS                            }....................
@type_check
def is_point(*points: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences are **two-dimensional points** (i.e.,
    contain exactly two numbers providing the X and Y coordinates of these
    points).

    Parameters
    ----------
    points : tuple[SequenceTypes]
        Tuple of all sequences to be tested.

    Returns
    ----------
    bool
        ``True`` only if these sequences are all two-dimensional points.
    '''

    # For efficiency, avoid testing whether these items are actually numbers.
    # While we could certainly do so, doing so is largely unhelpful; points
    # failing to contain numbers will generate readable exceptions elsewhere.
    return all(len(point) == 2 for point in points)


def is_left_of_vector(
    subject_point: SequenceTypes,
    vector_head_point: SequenceTypes,
    vector_tail_point: SequenceTypes,
) -> bool:
    '''
    ``True`` only if the passed two-dimensional subject point is spatially
    situated to the left of the two-dimensional vector defined by the passed
    two-dimensional head and tail points.

    Equivalently, this function returns ``True`` only if:

    * This vector is spatially situated to the right of this subject point.
    * A clockwise rotation rotates this subject point onto this vector.
    * A counter-clockwise rotation rotates this vector onto this subject point.
    * The determinant and hence sign of the cross product of this vector and
      this **subject vector** (i.e., vector whose head is this subject point and
      tail is the passed tail point) is positive.

    Derivation
    ----------
    We now derive the implementation of this function from the prior statement.
    If the angle between this vector and this subject vector (in that order) is:

    * Positive, a counter-clockwise rotation rotates the former onto the latter.
      In this case, this vector is spatially situated to the right of this
      subject point.
    * Negative, a clockwise rotation rotates the former onto the latter.  In
      this case, this vector is spatially situated to the left of this subject
      point.

    This function thus reduces to testing the sign of the angle between this
    vector and this subject vector. While there exist numerous means of doing
    so, the most computationally efficient in both space and time is to test the
    sign of the determinant of the cross product of these vectors -- requiring
    exactly four subtractions, two multiplications, and one comparison.

    The traditional definition of the cross product ``u × v`` of two arbitrary
    vectors ``u`` and ``v`` is given by the following two equivalent equalities:

    .. code::

       u × v = |u||v|sin(θ) n
               | i  j  k|
               |ux uy uz|
       u × v = |vx vy vz|

    Where:

    * ``|u|`` and ``|v|`` are the magnitudes of these vectors.
    * ``θ`` is the angle between these vectors.
    * ``n`` is a unit vector perpendicular to the plane containing these vectors
      in the direction oriented by the right-hand rule.
    * ``i``, ``j``, and ``k`` are the standard basis vectors.
    * ``ux``, ``uy``, and ``uz`` are the X, Y, and Z components of ``u``: e.g.,
      ``u = ux i + uy j + uz k``.
    * ``vx``, ``vy``, and ``vz`` are the X, Y, and Z components of ``v``: e.g.,
      ``v = vx i + vy j + vz k``.

    In this case, both ``u`` and ``v`` are two-dimensional, so:

    * ``n = k`` (i.e., ``u × v`` is a vector purely in the Z dimension).
    * ``uz = vz = 0``.

    By cofactor expansion, the latter equality for ``u × v`` then expands to:

    .. code::

               |uy uz|    |ux uz|    |ux uy|
       u × v = |vy vz|i - |vx vz|j + |vx vy|k
       u × v = (uy·vz - uz·vy)i - (ux·vz - uz·vx)j + (ux·vy - uy·vx)k
       u × v = (ux·vy - uy·vx)k

    By definition, the determinant ``det(u × v)`` of this cross product may be
    shown to similarly reduce to:

    .. code::

       det(u × v) = ux·vy - uy·vx

    The latter equality for ``u × v`` then reduces to:

    .. code::

       u × v = det(u × v) k

    The sign of the angle ``θ`` between ``u`` and ``v`` is then obtained by
    equating these two equalities and solving for ``θ`` as follows:

    .. code::

       |u||v|sin(θ) k = det(u × v) k
                det(u × v)
                ----------
       sin(θ) = |u||v|
                 [det(u × v)]
                 [----------]
       θ = arcsin[|u||v|    ]

    This computationally expensive equality may be reduced by noting that we
    only require the sign of ``θ``; the actual value of ``θ`` is irrelevant. In
    particular:

    * The ``arcsin()`` function is a monotonically increasing function vaguely
      resembling a wavy varient of the simple line ``y = x``, implying the sign
      of this function is *always* the sign of its operand, implying this
      function is ignorable with respect to this sign.
    * The quantity ``1/(|u||v|)`` is *always* positive, implying this quantity
      is also ignorable with respect to this sign.

    The sign of ``θ`` is thus the sign of ``det(u × v)``, the only remaining
    quantity with unknown sign in the above equality. This sign is positive if
    the following conditionality holds:

    .. code::

       ux·vy - uy·vx = det(u × v) > 0
       ux·vy > uy·vx

    Since this function accepts points rather than vectors, the former must be
    translated to the latter. Let:

    * ``h`` be the passed head point.
    * ``t`` be the passed tail point.
    * ``p`` be the passed subject point.
    * ``v`` be the corresponding subject vector.

    The X and Y components of ``u`` and ``v``  may then be defined in terms of
    these points as follows:

    * ``ux = hx - tx``.
    * ``uy = hy - ty``.
    * ``vx = px - tx``.
    * ``vy = py - ty``.

    The above conditionality then expands in terms of these points to:

    .. code::

       (hx - tx)·(py - ty) > (hy - ty)·(px - tx)

    The implementation of this function trivially follows.

    Parameters
    ----------
    vector_head_point : SequenceTypes
        2-sequence of the X and Y coordinates of the two-dimensional head point
        of this vector such that:
        * The first item is the X coordinate of this point.
        * The second item is the Y coordinate of this point.
    vector_tail_point : SequenceTypes
        2-sequence of the X and Y coordinates of the two-dimensional tail point
        of this vector.
    subject_point : SequenceTypes
        2-sequence of the X and Y coordinates of the two-dimensional **subject
        point** (i.e., point to test the orientation of this vector against).

    Returns
    ----------
    bool
        ``True`` only if this vector is to the right of this subject point.
    '''

    # If any passed sequence is *NOT* a point, raise an exception.
    die_unless_point(subject_point, vector_head_point, vector_tail_point)

    # See the above derivation for details. Math, you win all the efficiency.
    return (
        (vector_head_point[0] - vector_tail_point[0])*
            (subject_point[1] - vector_tail_point[1]) >
        (vector_head_point[1] - vector_tail_point[1])*
            (subject_point[0] - vector_tail_point[0])
    )

# ....................{ INTERSECTERS                       }....................
@type_check
def intersect_lines(
    line1_point1: SequenceTypes,
    line1_point2: SequenceTypes,
    line2_point1: SequenceTypes,
    line2_point2: SequenceTypes,
) -> SequenceTypes:
    '''
    Two-dimensional point intersecting the pair of lines defined by the passed
    pairs of two-dimensional points if any *or* raise an exception if these
    lines either infinitely overlap (i.e., are collinear) or never overlap
    (i.e., are parallel).

    Derivation
    ----------
    We now derive the implementation of this function from the two-point form of
    a two-dimensional line. Let:

    * ``c`` be the first passed point defining the first line.
    * ``d`` be the second passed point defining the first line.
    * ``s`` be the first passed point defining the second line.
    * ``t`` be the second passed point defining the second line.
    * ``e`` be any point residing on the first line.
    * ``u`` be any point residing on the second line.

    The intersection of these lines if any is the point satisfying the equality
    ``e = u``. For simplicity, let ``x = ex = ux`` and ``y = ex = ux`` be the
    X and Y coordinates of this point. The two-point form of the lines defined
    by these points constrained to intersect at this point is as follows:

    .. code::

                (dy - cy)(x - cx)
                -----------------
       y - cy =  dx - cx
                (ty - sy)(x - sx)
                ------------------
       y - sy =  tx - sx

    Determining this intersection point if any reduces to solving this system of
    two linear equations of two variables, which typically proceeds as follows:

    .. code::

       (y - cy)(dx - cx) = (dy - cy)(x - cx)
       (y - sy)(tx - sx) = (ty - sy)(x - sx)

       y·dx - y·cx - cy·dx + cy·cx = x·dy - dy·cx - x·cy + cy·cx
       y·tx - y·sx - sy·tx + sy·sx = x·ty - ty·sx - x·sy + sy·sx

       (cy - dy)x + (dx - cx)y = dx·cy - dy·cx
       (sy - ty)x + (tx - sx)y = tx·sy - ty·sx

    Let ``M`` be the coefficient matrix describing this linear system of
    equations such that:

    .. code::

           |cy-dy  dx-cx|
       M = |sy-ty  tx-sx|

    The determinant ``det(M)`` of this matrix is given by:

    .. code::

       det(M) = (cy-dy)·(tx-sx) - (dx-cx)·(sy-ty)

    This linear system of equations is then equivalent to this matrix equation:

    .. code::

        |x|   |dx·cy-dy·cx|
       M|y| = |tx·sy-ty·sx|

    This linear system of equations has a unique solution (implying these lines
    to intersect in a single point) if and only if:

    .. code::

       det(M) = (cy-dy)·(tx-sx) -  (dx-cx)·(sy-ty) != 0
                (cy-dy)·(tx-sx) != (dx-cx)·(sy-ty)

    Suppose ``det(M) != 0``, in which case this linear system of equastions has
    a unique solution. While there exist various (equally valid) means of
    obtaining this solution, the simplest is by invocation of Cramer's Rule.
    Let ``Mx`` be the matrix formed by replacing the first column of the
    coefficient matrix with the constant column vector and ``My`` the matrix
    formed by replacing the second column of the coefficient matrix with the
    constant column vector such that:

    .. code::

            |dx·cy-dy·cx  dx-cx      |
       Mx = |tx·sy-ty·sx  tx-sx      |
            |cy-dy        dx·cy-dy·cx|
       My = |sy-ty        tx·sy-ty·sx|

    The determinants ``det(Mx)`` and ``det(My)`` of these matrices are given by:

    .. code::

       det(Mx) = (dx·cy-dy·cx)(tx-sx) - (tx·sy-ty·sx)(dx-cx)
       det(My) = (tx·sy-ty·sx)(cy-dy) - (dx·cy-dy·cx)(sy-ty)

    Cramer's Rule then yields the X and Y coordinates of this point of
    intersection (with respect to these determinants) as follows:

    .. code::

           det(Mx)   (dx·cy-dy·cx)(tx-sx) - (tx·sy-ty·sx)(dx-cx)
           ------- = -------------------------------------------
       x = det(M)    (cy-dy)·(tx-sx) - (dx-cx)·(sy-ty)
           det(My)   (tx·sy-ty·sx)(cy-dy) - (dx·cy-dy·cx)(sy-ty)
           ------- = -------------------------------------------
       y = det(M)    (cy-dy)·(tx-sx) - (dx-cx)·(sy-ty)

    For efficiency, this solution is often rewritten into the equivalent form:

           (dx·cy-dy·cx)(tx-sx) - (tx·sy-ty·sx)(dx-cx)
           -------------------------------------------
       x = (dx-cx)·(ty-sy) - (dy-cy)·(tx-sx)
           (dx·cy-dy·cx)(ty-sy) - (tx·sy-ty·sx)(dy-cy)
           -------------------------------------------
       y = (dx-cx)·(ty-sy) - (dy-cy)·(tx-sx)

    The implementation of this function trivially follows.

    Parameters
    ----------
    line1_point1 : SequenceTypes
        2-sequence of the X and Y coordinates of the first two-dimensional point
        residing on the first of these lines (in any order) such that:
        * The first item is the X coordinate of this point.
        * The second item is the Y coordinate of this point.
    line1_point2 : SequenceTypes
        2-sequence of the X and Y coordinates of the second two-dimensional
        point residing on the first of these lines (in any order).
    line2_point1 : SequenceTypes
        2-sequence of the X and Y coordinates of the first two-dimensional
        point residing on the second of these lines (in any order).
    line2_point2 : SequenceTypes
        2-sequence of the X and Y coordinates of the second two-dimensional
        point residing on the second of these lines (in any order).

    Returns
    ----------
    SequenceTypes
        2-sequence of the X and Y coordinates of the two-dimensional point
        intersecting the pair of lines defined by these two-dimensional points.
        The type of this sequence is the same as that of the first passed point.

    Raises
    ----------
    BetseMathLineException
        If these lines either:
        * Infinitely overlap (i.e., are collinear).
        * Never overlap (i.e., are parallel).
    '''

    # If any passed sequence is *NOT* a point, raise an exception.
    die_unless_point(
        line1_point1, line1_point2, line2_point1, line2_point2)

    # X coordinates of these points in the nomenclature documented above.
    cx = line1_point1[0]
    dx = line1_point2[0]
    sx = line2_point1[0]
    tx = line2_point2[0]

    # Y coordinates of these points in the nomenclature documented above.
    cy = line1_point1[1]
    dy = line1_point2[1]
    sy = line2_point1[1]
    ty = line2_point2[1]

    # Run and rise of the slope of the first line.
    dcx = dx - cx
    dcy = dy - cy

    # Run and rise of the slope of the second line.
    tsx = tx - sx
    tsy = ty - sy

    # Determinant of the coefficient matrix solving this intersection.
    detM = dcx*tsy - dcy*tsx

    # If this determinant is zero, no unique solution exists, in which case
    # these lines either infinitely or never intersect. Since this constitutes
    # an error condition in either case, raise an exception. While we could
    # detect and differentiate these two error conditions in this exception
    # message, a lazy approach currently bests an accurate approach.
    if detM == 0:
        raise BetseMathLineException(
            'Line through points {!r} and {!r} and '
            'line through points {!r} and {!r} '
            'either infinitely or never intersect.'.format(
                cx, dx, sx, tx))
    # Else, this determinant is non-zero, in which a unique solution exists and
    # these lines uniquely intersect at a single point, which we now obtain.

    #                                  |dx dy|
    # Determinant of the ad-hoc matrix |cx cy|.
    detdc = dx*cy - dy*cx

    #                                  |tx ty|
    # Determinant of the ad-hoc matrix |sx sy|.
    detts = tx*sy - ty*sx

    # 2-tuple of the X and Y coordinates of this intersection, calculated via
    # Cramer's Rule from the requisite determinants.
    intersection_point = (
        (detdc*tsx - detts*dcx) / detM,
        (detdc*tsy - detts*dcy) / detM,
    )

    # Return this intersection as a sequence of the same type as the first
    # passed point.
    return iterables.to_iterable(
        iterable=intersection_point, cls=type(line1_point1))
