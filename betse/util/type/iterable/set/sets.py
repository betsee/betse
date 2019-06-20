#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **set** (i.e., :class:`set`-like types or instances) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import type_check, SetType

# ....................{ METACLASSES                       }....................
@type_check
def symmetric_difference(*sets: SetType) -> SetType:
    '''
    **Symmetric difference** (i.e., set of all objects in only a single passed
    set) of all passed sets.

    This function generalizes the :meth:`set.symmetric_difference` method to an
    arbitrary number of sets. Whereas the :meth:`set.intersection` and
    :meth:`set.union` methods natively support an arbitrary number of sets, the
    :meth:`set.symmetric_difference` method does *not*. This function corrects
    this unfortunate oversight.

    See Also
    ----------
    https://bugs.python.org/issue17854
        Official Python discussion of this issue.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import sequences

    #FIXME: Raise an exception if less than two sets were passed!
    # If no sets were passed, raise an exception.
    sequences.die_if_empty(sets, label='Set')

    # Symmetric difference of these sets to be returned, initialized to a
    # shallow copy of the first such set.
    set_difference = sets[0].copy()

    # For each subsequent such set...
    for set_next in sets[1:]:
        # Reduce this symmetric difference against this set.
        set_difference = set_difference.symmetric_difference(set_next)

    # Return this symmetric difference.
    return set_difference
