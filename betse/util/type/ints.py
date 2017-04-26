#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level integer facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseIntException
from betse.util.type import types
from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
BITS_PER_BYTE = 8
'''
Number of bits per byte.
'''

# ....................{ CONSTANTS ~ max                    }....................
BYTE_VALUE_MAX = 255
'''
Maximum value of unsigned bytes.
'''


INT_VALUE_MAX_32_BIT = 1 << 32
'''
Maximum value for integer variables of internal type `Py_ssize_t` on 32-bit
Python interpreters.

This value is suitable for comparison with `sys.maxsize`, the maximum value of
these variables on the current system.
'''

# ....................{ CONSTANTS ~ units ~ size           }....................
# Size denominations in base 2 rather than base 10, for mild efficiency.

KiB = 1 << 10
'''
Number of bytes in one kibibyte (i.e., `2**10`).
'''


MiB = 1 << 20
'''
Number of bytes in one mebibyte (i.e., `2**20`).
'''


GiB = 1 << 30
'''
Number of bytes in one gibibyte (i.e., `2**30`).
'''


TiB = 1 << 40
'''
Number of bytes in one tebibyte (i.e., `2**40`).
'''

# ....................{ CONSTANTS ~ units ~ inverse        }....................
INVERSE_CENTI = 1e2
'''
Inverse of the centi- unit prefix (i.e., `10**−2`), typically used as a
multiplicative factor for upscaling quantities from centi-prefixed units to
unprefixed units.
'''


INVERSE_MILLI = 1e3
'''
Inverse of the milli- unit prefix (i.e., `10**−3`), typically used as a
multiplicative factor for upscaling quantities from milli-prefixed units to
unprefixed units.
'''


INVERSE_MICRO = 1e6
'''
Inverse of the micro- unit prefix (i.e., `10**−6`), typically used as a
multiplicative factor for upscaling quantities from micro-prefixed units to
unprefixed units.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless(*objects) -> None:
    '''
    Raise an exception unless all passed objects are integers.

    Parameters
    ----------
    objects : tuple
        Tuple of all objects to be validated.

    Raises
    ----------
    BetseIntException
        If any passed object is _not_ an integer.
    '''

    for obj in objects:
        if not types.is_int(obj):
            raise BetseIntException(
                'Object "{}" not an integer.'.format(obj))


@type_check
def die_unless_positive(*numbers: int, label: str = 'Integer') -> None:
    '''
    Raise an exception prefixed by the passed label unless all passed objects
    are positive integers.

    Parameters
    ----------
    numbers : tuple
        Tuple of all integers to be validated.
    label : optional[str]
        Human-readable label prefixing exception messages raised by this method.
        Defaults to a general-purpose string.

    Raises
    ----------
    BetseIntException
        If any passed object is _not_ a positive integer.
    '''

    # For each passed integer...
    for number in numbers:
        # If this integer is non-positive, raise an exception.
        if number <= 0:
            raise BetseIntException(
                '{} "{}" not positive.'.format(label.capitalize(), number))

# ....................{ TESTERS                            }....................
@type_check
def is_even(number: int) -> bool:
    '''
    ``True`` only if the passed integer is even.
    '''

    return number % 2 == 0


@type_check
def is_odd(number: int) -> bool:
    '''
    ``True`` only if the passed integer is odd.
    '''

    return number % 2 == 1
