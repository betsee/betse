#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **bitwise** (i.e., callables operating on one or more individual bits
of integer-based bit fields) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import type_check

# ....................{ TESTERS                           }....................
@type_check
def is_bit_on(bit_field: int, bit_mask: int) -> bool:
    '''
    ``True`` only if the bit in the passed bit field uniquely identified by the
    passed bit mask is **on** (i.e., one rather than zero).

    Parameters
    ----------
    bit_field : int
        **Bit field** (i.e., integer whose individual bits collectively
        comprise a rudimentary array of boolean flags) to be tested.
    bit_mask : int
        **Bit mask** (i.e., integer containing exactly one non-zero bit at the
        position uniquely identifying the desired bit) of the bit to be tested.

    Returns
    ----------
    bool
        ``True`` only if this bit in this bit field is on.
    '''

    # Return true only if this bit in this bit field is on. Dismantled, this
    # is:
    #
    # * "bit_field & bit_mask", either:
    #   * If this bit is on in this bit field, "bit_mask" and hence non-zero.
    #   * If this bit is off in this bit field, zero.
    # * "bool(...)", reducing this integer to a boolean.
    return bool(bit_field & bit_mask)
