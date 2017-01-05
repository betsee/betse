#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **numeric** (i.e., number in the general sense and hence equally
applicable to both integers and floats) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseNumericException
from betse.util.type.types import type_check, IterableTypes, NumericTypes

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

# ....................{ SUMMERS                            }....................
@type_check
def sum_by_index(iterable: IterableTypes, subiterable_index: object) -> (
    NumericTypes):
    '''
    Summation of each **number** (i.e., integer or float) at the passed key or
    index of each subiterable of the passed iterable.

    Each number at the passed key or index of each subiterable of this iterable
    is summed with the `+` operator, implicitly calling the `__sum__()`
    special method of these numbers. The type of the resulting number is the
    widest type of these numbers. Specifically:

    * If any such number is a float, the returned number is also a float.
    * Else, the returned number is an integer.

    For disambiguity, each number is ideally but _not_ necessarily of the same
    numeric type: namely, either all integers or all floats.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of subiterables to be summed.
    subiterable_index : object
        Object with which to index each subiterable of this iterable. The type
        of this object _must_ be a type accepted by the `__getitem__()` special
        method of each subiterable. Specifically, if each subiterable is a:
        * **Mapping** (e.g., :class:`dict`), this object _must_ be hashable.
        * **Sequence** (e.g., :class:`list`, :class:`tuple`), this object
          _must_ be either:
          * An integer.
          * A :func:`slice` object.

    Returns
    ----------
    NumericTypes
        Summation of each number at this key or index of each subiterable of
        this iterable.
    '''

    # Efficiency and simplicity combine here to form MegaFastSimple.
    return sum(subiterable[subiterable_index] for subiterable in iterable)
