#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Primitive vector functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, SequenceTypes

# ....................{ TESTERS                            }....................
#FIXME: Implement us up.
#FIXME: Document the derivation.
#FIXME: Unit test us up.
@type_check
def is_right_of_point(
    vector_head_point: SequenceTypes,
    vector_tail_point: SequenceTypes,
    subject_point: SequenceTypes,
) -> bool:
    '''
    ``True`` only if the vector defined by the passed head and tail points is
    spatially situated to the right of the passed subject point.

    Equivalently, this function returns ``True`` only if the:

    * This subject point is spatially situated to the left of this vector.
    * The cross product of this vector and the vector whose head is this subject
      point and whose tail is the passed tail point is positive.

    Derivation
    ----------

    Parameters
    ----------
    vector_head_point : SequenceTypes
        2-sequence of the X and Y coordinates of the head point of this vector
        such that:
        * The first item is the X coordinate of this point.
        * The second item is the Y coordinate of this point.
    vector_tail_point : SequenceTypes
        2-sequence of the X and Y coordinates of the tail point of this vector.
    subject_point : SequenceTypes
        2-sequence of the X and Y coordinates of the **subject point** (i.e.,
        point to test the orientation of this vector against).

    Returns
    ----------
    bool
        ``True`` only if this vector is to the right of this subject point.
    '''

    pass
