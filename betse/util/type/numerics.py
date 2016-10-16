#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **numeric** (i.e., number in the general sense and hence equally
applicable to both integers and floats) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseNumericException
from betse.util.type.types import type_check, NumericTypes

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_positive(*numbers: NumericTypes, label: str = 'Number') -> None:
    '''
    Raise an exception prefixed by the passed label unless all passed objects
    are positive integers or floats.

    Parameters
    ----------
    numbers : tuple
        Tuple of all objects to be validated.
    label : optional[str]
        Human-readable label prefixing exception messages raised by this method.
        Defaults to a general-purpose string.

    Raises
    ----------
    BetseNumericException
        If any passed object is neither a positive integer nor float.
    '''

    # For each passed number...
    for number in numbers:
        # If this number is non-positive...
        if number <= 0:
            # Raise an exception.
            raise BetseNumericException(
                '{} "{}" not positive.'.format(label.capitalize(), number))
