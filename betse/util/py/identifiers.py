#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **Python identifier** (i.e., class, module, or attribute name)
facilities.
'''

# ....................{ IMPORTS                            }....................
import re
from betse.exceptions import BetseStrException
from betse.util.type.types import type_check, RegexMatchType

# ....................{ CLASSES                            }....................
IDENTIFIER_ALPHANUMERIC_CHAR_CLASS = r'a-zA-Z0-9'
'''
Character class (excluding `[` and `]` delimiters) matching any character of a
**Python identifier** (i.e., class, module, or attribute name), excluding the
underscore character.

Of necessity, this character class is equivalent to the character class of all
alphanumeric ASCII characters.
'''


IDENTIFIER_CHAR_CLASS = IDENTIFIER_ALPHANUMERIC_CHAR_CLASS + r'_'
'''
Character class (excluding `[` and `]` delimiters) matching any character of a
**Python identifier** (i.e., class, module, or attribute name).
'''

# ....................{ REGEXES                            }....................
IDENTIFIER_CHAR_CAMELCASE_REGEX = (
    r'(?:(?<=[a-z0-9])|(?!^))([A-Z])(?=[a-z]|$)')
'''
Uncompiled regular expression matching the next first character of a contiguous
run of uppercase characters in a CamelCase-formatted Python identifier (e.g.,
the `I` in `Capitel IV`), excluding the first character of such identifier..

Examples
----------

This expression is intended to be used in substitutions converting CamelCase to
some other format. For example, to convert CamelCase to snake_case:

    >>> from betse.util.type import regexes
    >>> regexes.replace_substrs(
    ...     'MesseIoXaVIaX',
    ...     rexeges.IDENTIFIER_CHAR_CAMELCASE_REGEX,
    ...     r'_\1')
    Messe_io_xa_via_x
'''


IDENTIFIER_UNQUALIFIED_REGEX = r'[{}]+'.format(IDENTIFIER_CHAR_CLASS)
'''
Uncompiled regular expression matching an **unqualified Python identifier**
(i.e., module or attribute name _not_ prefixed by a package or module name).
'''


IDENTIFIER_QUALIFIED_REGEX = (
    r'(?:{identifier_unqualified}\.)*{identifier_unqualified}'.format(
        identifier_unqualified=IDENTIFIER_UNQUALIFIED_REGEX))
'''
Uncompiled regular expression matching a **qualified Python identifier** (i.e.,
module or attribute name possibly prefixed by a package or module name).
'''

# ....................{ PRIVATE ~ regexes                  }....................
# This regular expression is frequently referenced at application startup and
# thus unconditionally compiled. Specifically, this regex is referenced by the
# sanitize_camelcase() function called by the var_alias() descriptor
# repeatedly instantiated at class scope by "Parameters"-related classes.
_IDENTIFIER_SANITIZE_CAMELCASE_REGEX = re.compile(
    r'(?:^|[^a-zA-Z0-9]+)([a-z])?')
'''
Compiled regular expression matching either the string start *or* one or more
alphanumeric ASCII characters optionally followed by a lowercase alphabetic
ASCII character, internally required by the :func:`sanitize_camelcase` function.
'''

# ....................{ CONVERTERS ~ camel                 }....................
def convert_camelcase_to_snakecase(identifier: str) -> str:
    '''
    Convert the passed CamelCase-formatted Python identifier to snake_case
    (e.g., from ``ThePMRC`` to ``the_pmrc``).
    '''

    # Avoid circular import dependencies.
    from betse.util.type import regexes

    # Munge this identifier.
    return regexes.replace_substrs(
        text=identifier,
        regex=IDENTIFIER_CHAR_CAMELCASE_REGEX,
        replacement=r'_\1',
    ).lower()


def convert_camelcase_to_whitespaced_lowercase(identifier: str) -> str:
    '''
    Convert the passed CamelCase-formatted Python identifier to whitespaced
    lowercase (e.g., from `CleanseIII` to `cleanse iii`).
    '''

    # Avoid circular import dependencies.
    from betse.util.type import regexes

    # Munge this identifier.
    return regexes.replace_substrs(
        text=identifier,
        regex=IDENTIFIER_CHAR_CAMELCASE_REGEX,
        replacement=r' \1',
    ).lower()

# ....................{ SANITIZERS ~ camelcase             }....................
@type_check
def sanitize_camelcase(identifier: str) -> str:
    '''
    Sanitize the passed string into a valid Python identifier in CamelCase
    format derived from this string (e.g., from ``sim-gnr`` to ``SimGnr``).

    Specifically, this function replaces each substring of one or more
    alphanumeric ASCII characters with:

    * If this substring is followed by a lowercase alphabetic ASCII
      character, that character uppercased (e.g., from ``['spider jerusalem']``
      to ``SpiderJerusalem``).
    * Else, the empty string (e.g., from ``['']`` to ````).

    Parameters
    ----------
    identifier : str
        String to be converted into a valid CamelCase-style Python identifier.

    Returns
    ----------
    str
        Valid CamelCase-style Python identifier converted from this string.

    Raises
    ----------
    BetseStrException
        If the passed string contains no alphanumeric characters, in which case
        the resulting identifier would be empty and hence invalid.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import regexes

    # Identifier sanitized from this unsanitized identifier as detailed above.
    identifier_sanitized = regexes.replace_substrs(
        text=identifier,
        regex=_IDENTIFIER_SANITIZE_CAMELCASE_REGEX,
        replacement=_sanitize_camelcase_match,
    )

    # Return this sanitized identifier is empty, this unsanitized identifier
    # must necessarily have contained no alphanumeric characters. In this case,
    # raise an exception.
    if not identifier_sanitized:
        raise BetseStrException(
            'Unsanitized identifier "{}" contains no '
            'alphanumeric characters.'.format(identifier))

    # Return this sanitized identifier.
    return identifier_sanitized


def _sanitize_camelcase_match(match: RegexMatchType) -> str:
    '''
    Uppercase the lowercase alphabetic ASCII character optionally grouped by the
    passed regular expression match object if any *or* return the empty string
    otherwise.

    Parameters
    ----------
    match: RegexMatchType
        Match object optionally grouping a lowercase alphabetic ASCII character.

    Returns
    ----------
    str
        Either:
        * If this match object grouped a lowercase alphabetic ASCII character,
          this character uppercased.
        * Else, the empty string.
    '''

    return (match.group(1) or '').upper()

# ....................{ SANITIZERS ~ snakecase             }....................
#FIXME: Generalize to sanitize arbitrary strings. See the general-purpose
#sanitize_camelcase() function for inspiration.
@type_check
def sanitize_snakecase(identifier: str) -> str:
    '''
    Sanitize the passed string into a valid Python identifier in snake_case
    format derived from this string (e.g., from ``sim-gnr`` to ``sim_gnr``).

    Specifically, this function:

    * Replaces each of the following characters in this string with an
      underscore:
      * A hyphen.

    Parameters
    ----------
    identifier : str
        String to be converted into a valid snake_case-style Python identifier.

    Returns
    ----------
    str
        Valid snake_case-style Python identifier converted from this string.
    '''

    return identifier.replace('-', '_')
