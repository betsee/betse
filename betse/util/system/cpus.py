#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level Central Processing Unit (CPU) facilities.

Caveats
----------
CPU-specific functions (e.g., `is_64_bit()`) are considered poor form. Call such
functions only when absolutely necessary.
'''

# ....................{ IMPORTS                            }....................
import sys

# ....................{ TESTERS ~ arch                     }....................
def is_32_bit():
    '''
    True if the current CPU architecture is 32-bit.
    '''
    return not is_64_bit()

def is_64_bit():
    '''
    True if the current CPU architecture is 64-bit.
    '''
    # Avoid circular import dependencies.
    from betse.util.type import ints

    # There exist several alternative means of testing the same condition: e.g.,
    #
    #     return 'PROCESSOR_ARCHITEW6432' in os.environ
    #
    # The current approach, however, is more portable and hence ideal.
    return sys.maxsize > ints.INT_VALUE_MAX_32_BIT

# --------------------( WASTELANDS                         )--------------------
#FUXME: Research us up. Is this *REALLY* the canonical means of determining
#this under Python3? (This seems fairly terrible, honestly. Is this actually
#cross-platform-portable?)
