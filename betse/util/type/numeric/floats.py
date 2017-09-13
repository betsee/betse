#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level floating point facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.call.memoizers import callable_cached
from betse.util.type.types import type_check, RegexCompiledType

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
        >>> floats.get_precision(0.110001000000000009)
        17
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # String formatted from this number, guaranteed by Python internals to
    # comply with one of the following two formats:
    #
    # * Decimal notation (e.g., "3.1415").
    # * Scientific notation (e.g., "3.1415e+92", "3.141592e-65").
    number_str = str(number)

    # Compiled regex matching a floating point number in scientific notation.
    scientific_notation_regex = _get_scientific_notation_regex()

    # If this number is in scientific notation
    scientific_notation_groups = regexes.get_match_groups_numbered_if_any(
        text=number_str, regex=scientific_notation_regex)

    # If this number is in decimal notation...
    if scientific_notation_groups is None:
        # 0-based index of the decimal in this number's string.
        decimal_index = number_str.index('.')

        # Substring of this number's fractional part excluding this decimal.
        decimals = number_str[decimal_index + 1:]

        # This number's precision is the number of digits in this substring.
        precision = len(decimals)
    # Else, this number is in scientific notation. In this case..
    else:
        # This number's precision is the absolute value (i.e., ignoring the
        # sign) of this number's exponent in this notation.
        precision = abs(int(scientific_notation_groups[1]))

    # Return this precision.
    return precision

# ....................{ GETTERS ~ private                  }....................
@callable_cached
def _get_scientific_notation_regex() -> RegexCompiledType:
    '''
    Compiled regular expression matching a floating point number in scientific
    notation encapsulated as a string (e.g., ``6.66e9``), capturing the
    significand and exponent of this number into the first two match groups.

    For efficiency in downstream clients (e.g., BETSEE) frequently calling the
    :func:`is_ml` function and hence requiring this expression be compiled, this
    expression is intentionally pre-compiled rather than uncompiled and thus
    returned as a cached getter.

    See Also
    ----------
    https://stackoverflow.com/a/658662/2809027
        StackOverflow answer strongly inspiring this regular expression.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Create, return, and cache this expression.
    return regexes.compile(
        # Significand (captured).
        r'('
            # Negative prefix (optional).
            r'-?'
            # Non-fractional part of the significand. Either:
            r'(?:'
                # Zero.
                r'0|'
                # One or more digits whose first digit is *NOT* zero.
                r'[1-9]\d*'
            r')'
            # Fractional part of the significand (optional).
            r'(?:\.\d*)?'
        r')'
        # Exponent prefix.
        r'[eE]'
        # Exponent (captured).
        r'('
            # Exponent sign.
            r'[+\-]'
            # Exponent digits.
            r'\d+'
        r')'
    )
