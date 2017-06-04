#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive two-dimensional point functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMathException
from betse.util.type.types import type_check, SequenceTypes

# ....................{ EXCEPTIONS                         }....................
def die_unless_point(*points: SequenceTypes) -> None:
    '''
    Raise an exception unless all passed sequences contain exactly two numbers.

    Parameters
    ----------
    points : tuple[SequenceTypes]
        Tuple of all sequences to be validated.

    Raises
    ----------
    BetseMathException
        If any such sequence does *not* contain exactly two numbers.
    '''

    if not is_point(*points):
        for point in points:
            if not is_point(point):
                raise BetseMathException(
                    'Sequence not a two-dimensional point '
                    '(i.e., length != 2): {!r}'.format(point))

# ....................{ TESTERS                            }....................
@type_check
def is_point(*points: SequenceTypes) -> bool:
    '''
    ``True`` only if all passed sequences contain exactly two numbers,
    presumably signifying the X and Y coordinates (in order) of points.

    Parameters
    ----------
    points : tuple[SequenceTypes]
        Tuple of all sequences to be tested.

    Returns
    ----------
    bool
        ``True`` only if these sequences all contain exactly two numbers.
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

    However, this function accepts points rather than vectors, necessitating
    translation from the former to the latter. Denote ``h`` the passed head
    point, ``t`` the passed tail point, ``p`` the passed subject point, and
    ``v`` the corresponding subject vector. Given these points, the vector
    components of ``u`` and ``v``  are defined as follows:

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
