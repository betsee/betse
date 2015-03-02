#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level regular expression (regex) facilities.
'''

# ....................{ IMPORTS                            }....................
from sre_parse import Pattern
import re

# ....................{ CONSTANTS ~ python                 }....................
PYTHON_IDENTIFIER_CHAR_CLASS = r'a-zA-Z0-9_'
'''
Character class (excluding `[` and `]` delimiters) matching any character of a
**Python identifier** (i.e., class, function, module, or variable name).
'''

PYTHON_IDENTIFIER_REGEX_RAW = r'[{}]+'.format(PYTHON_IDENTIFIER_CHAR_CLASS)
'''
Uncompiled regular expression matching a **Python identifier** (i.e., class,
function, module, or variable name).
'''

# ....................{ REPLACERS                          }....................
def remove_substrings(text: str, regex, **kwargs) -> str:
    '''
    Remove all substrings in the passed string matching the passed regular
    expression.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object).

    This function accepts the same optional keyword arguments as `re.sub()`.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        For further details on such regular expressions and keyword arguments.
    '''
    return substitute_substrings(text, regex, '', **kwargs)

# ....................{ SUBSTITUTERS                       }....................
def substitute_substrings(text: str, regex, substitution, **kwargs) -> str:
    '''
    Substitute all substrings in the passed string matching the passed regular
    expression with the passed substitution.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object). Likewise, such substitution may be
    either a string *or* function.

    This function accepts the same optional keyword arguments as `re.sub()`.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        For further details on such regular expressions, substitutions, and
        keyword arguments.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    return re.sub(regex, substitution, text, **kwargs)

# --------------------( WASTELANDS                         )--------------------
    # if isinstance(regex, Pattern):
    # assert isinstance(substitution, str), '"{}" not a string.'.format(
    #     substitution)
