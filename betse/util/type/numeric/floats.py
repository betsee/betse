#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level floating point facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.type.call.memoizers import callable_cached
from betse.util.type.types import type_check, RegexCompiledType

# ....................{ TESTERS                            }....................
@type_check
def is_float_str(text: str) -> bool:
    '''
    ``True`` only if the passed string syntactically conforms to either the
    decimal format (e.g., ``6.669``) or scientific notation format (e.g.,
    ``6.66e9``) expected by floating point numbers.

    Equivalently, this tester returns ``True`` only if this string is losslessly
    convertable into a floating point number.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Return True only if this string matches a floating point format.
    return regexes.is_match(text=text, regex=get_float_regex())

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
    number_text = str(number)

    # Scientific notation match groups if this number is in scientific notation,
    # capturing the significand and exponent of this number into these groups.
    number_text_groups = regexes.get_match_groups_named(
        text=number_text, regex=get_float_regex())

    # If this number is in decimal notation...
    if number_text_groups['exponent'] is None:
        # Fractional part of this number if any *OR* the empty string otherwise.
        # Retrieving this substring is complicated by the fact that one of two
        # alternate match groups may yield this substring at a time.
        significand_frac = (
            number_text_groups['significand_frac_empty'] or
            number_text_groups['significand_frac_nonempty'] or '')
        logs.log_debug('significand (fractional): %s', significand_frac)
        # print('significand (fractional): %s' % significand_frac)

        # This number's precision is the number of digits in this part.
        precision = len(significand_frac)
    # Else, this number is in scientific notation. In this case..
    else:
        # print('exponent: %s' % number_text_groups['exponent'])
        # This number's precision is the absolute value (i.e., ignoring the
        # sign) of this number's exponent in this notation.
        precision = abs(int(number_text_groups['exponent']))

    # Return this precision.
    return precision

# ....................{ GETTERS ~ regex                    }....................
# For efficiency in downstream clients (e.g., BETSEE) frequently calling the
# is_ml() function and hence requiring this expression be compiled, this
# expression is intentionally pre-compiled rather than uncompiled and thus
# returned as a cached getter.

@callable_cached
def get_float_regex() -> RegexCompiledType:
    '''
    Compiled regular expression matching a floating point number in either
    decimal notation (e.g., ``6.69``) *or* scientific notation (e.g., ``6.6e9``)
    encapsulated as a string.

    This expression captures the following named match groups (in order):

    1. ``significand``, yielding the mandatory significand of this number
       *including* optional sign prefix (e.g., ``-6.7`` given ``-6.7e-9``).
    1. ``significand_nonfrac_empty``, yielding the optional non-fractional
       digits of this significand *including* optional sign prefix (e.g., ``-6``
       given ``-6.7e-9``). This group's value is ``None`` for numbers prefixed
       by ``.`` rather than a digit (e.g., ``-.67e-9``).
    1. ``significand_frac_nonempty``, alternatively yielding the optional
       fractional digits of this significand in a manner guaranteed to be
       non-empty (e.g., ``None`` given ``-6.7e-9``).
    1. ``significand_frac_empty``, alternatively yielding the optional
       fractional digits of this significand in a manner *not* guaranteed to be
       non-empty (e.g., ``7`` given ``-6.7e-9``). This edge case is required to
       support floating point numbers of both the form
       ``{significand_nonfrac}.{significand_frac_empty}`` *and*
       ``.{significand_frac_nonempty}``. Since the standard :mod:`re` module
       does *not* support the ``(?|...)``-style branch reset syntax supported by
       the third-party :mod:`regex` module, this is the best we can do without
       burdening the codebase with *yet another* third-party dependency.
    1. ``exponent``, yielding the optional exponent of this number *including*
       optional sign prefix (e.g., ``-9`` given ``6.7e-9``).

    See Also
    ----------
    https://jdreaver.com/posts/2014-07-28-scientific-notation-spin-box-pyside.html
        Blog article partially inspiring this implementation.
    https://stackoverflow.com/a/658662/2809027
        StackOverflow answer partially inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Create, return, and cache this expression.
    return regexes.compile(
        # Significand (captured).
        r'(?P<significand>'
            # Negative or positive prefix (optional).
            r'[+-]?'
            # Either:
            r'(?:'
            # * Significand prefixed by one or more non-fractional digits.
                # Non-fractional digits.
                r'(?P<significand_nonfrac_empty>\d+)'
                # Fractional digits (optional).
                r'(?:\.(?P<significand_frac_empty>\d*))?'
            r'|'
            # * Significand containing only fractional digits.
                r'\.(?P<significand_frac_nonempty>\d+)'
            r')'
        r')'
        # Exponent (optional).
        r'(?:'
            # Exponent prefix.
            r'[eE]'
            # Exponent (captured).
            r'(?P<exponent>'
                # Exponent sign (optional).
                r'[+-]?'
                # Exponent digits.
                r'\d+'
            r')'
        r')?'
    )


#FIXME: Obsoleted by the above regular expression.
@callable_cached
def get_float_scientific_notation_regex() -> RegexCompiledType:
    '''
    Compiled regular expression matching a floating point number in scientific
    notation encapsulated as a string (e.g., ``6.6e9``).

    This expression captures the following positional match groups (in order):

    1. Significand of this number (e.g., ``6.6`` in the case of ``6.6e9``).
    1. Exponent of this number (e.g., ``9`` in the case of ``6.6e9``).

    See Also
    ----------
    https://stackoverflow.com/a/658662/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Create, return, and cache this expression.
    return regexes.compile(
        # Significand (captured).
        r'('
            # Significand sign (optional).
            r'[+-]?'
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
            # Exponent sign (optional).
            r'[+-]?'
            # Exponent digits.
            r'\d+'
        r')'
    )
