#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **character** (i.e., string of length 1) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseCharException
from betse.util.type.types import type_check

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_char(text: str) -> None:
    '''
    Raise an exception unless the passed string is a **character** (i.e., of
    length 1).

    Parameters
    ----------
    text : str
        String to be validated.

    Raises
    ----------
    BetseCharException
        If this string is either:
        * The empty string.
        * Of length greater than 1.
    '''

    # If this string is *NOT* a character, raise an exception.
    if not is_char(text):
        raise BetseCharException(
            'String "{}" not a character '
            '(i.e., either empty or of length >= 2).'.format(text))

# ....................{ TESTERS ~ case                     }....................
@type_check
def is_char(text: str) -> bool:
    '''
    ``True`` only if the passed string is a **character** (i.e., of length 1).

    Parameters
    ----------
    text : str
        String to be tested.

    Returns
    ----------
    bool
        ``True`` only if this string is a character.
    '''

    # Don't judge us.
    return len(text) == 1
