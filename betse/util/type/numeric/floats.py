#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level floating point facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check

# ....................{ GETTERS                            }....................
@type_check
def get_precision(number: float) -> int:
    '''
    Precision of the passed floating point number.

    Precision is defined as the length of this number's significand (excluding
    leading hidden digit 1), equivalent to the number of base-10 digits in the
    fractional segment of this number *after* accounting for approximation
    errors in floating point arithmetic.

    Examples
    ----------
        >>> from betse.util.type.numeric import floats
        >>> floats.get_precision(len(str(0.110001000000000009)) - 2
        17
    '''

    return len(str(number)) - 2
