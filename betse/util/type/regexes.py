#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level regular expression (regex) facilities.
'''

# ....................{ IMPORTS                            }....................
import re

# ....................{ CONSTANTS ~ python                 }....................
PYTHON_IDENTIFIER_CHAR_CLASS = r'a-zA-Z0-9_'
'''
Character class (excluding `[` and `]` delimiters) matching any character of a
**Python identifier** (i.e., class, function, module, or variable name).
'''

PYTHON_IDENTIFIER_UNQUALIFIED_REGEX_RAW = r'[{}]+'.format(
    PYTHON_IDENTIFIER_CHAR_CLASS)
'''
Uncompiled regular expression matching an **unqualified Python identifier**
(i.e., class, function, module, or variable name *not* prefixed by a package or
module name).
'''

PYTHON_IDENTIFIER_QUALIFIED_REGEX_RAW =\
    r'(?:{identifier_unqualified}\.)*{identifier_unqualified}'.format(
        identifier_unqualified = PYTHON_IDENTIFIER_UNQUALIFIED_REGEX_RAW)
'''
Uncompiled regular expression matching an **qualified Python identifier**
(i.e., class, function, module, or variable name possibly prefixed by a package
or module name).
'''

# ....................{ MATCHERS                           }....................
def get_match_groups_named(text: str, regex, **kwargs) -> list:
    '''
    Get the dictionary from explicitly named groups to substrings matched from
    the passed string with the passed regular expression.

    Unmatched groups will have the value `None`. Unnamed (i.e., strictly
    numeric) groups will *not* be provided by the returned dictionary,
    regardless of whether such groups matched. If this is undesirable, consider
    calling `get_match_groups_numeric()` instead.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object). This function accepts the same optional
    keyword arguments as `re.match()`.

    Returns
    ----------
    dict
        Dictionary mapping such groups.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''
    return get_match(text, regex, **kwargs).groupdict()

def get_match_groups_numeric(text: str, regex, **kwargs) -> list:
    '''
    Get the list of all groups matched from the passed string with the passed
    regular expression, ordered by the left-to-right lexical position at which
    each such group was matched.

    Unmatched groups will have the value `None`.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object). This function accepts the same optional
    keyword arguments as `re.match()`.

    Returns
    ----------
    list
        List of such groups.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''
    return get_match(text, regex, **kwargs).groups()

def get_match(text: str, regex, **kwargs):
    '''
    Get the match object obtained by matching the passed string with the passed
    regular expression.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object). This function accepts the same optional
    keyword arguments as `re.match()`.

    Returns
    ----------
    SRE_Match
        Such match object.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    return re.match(regex, text, **kwargs)

# ....................{ REPLACERS                          }....................
def remove_substrings(text: str, regex, **kwargs) -> str:
    '''
    Remove all substrings in the passed string matching the passed regular
    expression.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object). This function accepts the same optional
    keyword arguments as `re.sub()`.

    Returns
    ----------
    str
        Passed string containing no substrings matching such regular expression.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''
    return substitute_substrings(text, regex, '', **kwargs)

# ....................{ SUBSTITUTERS                       }....................
def substitute_substrings(text: str, regex, substitution, **kwargs) -> str:
    '''
    Substitute all substrings in the passed string matching the passed regular
    expression with the passed substitution.

    Such regular expression may be either a string *or* `Pattern` (i.e.,
    compiled regular expression object), while such substitution may be either a
    string *or* function. This function accepts the same optional keyword
    arguments as `re.sub()`.

    Returns
    ----------
    str
        Passed string with all substrings matching such regular expression
        replaced with such substitution.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    return re.sub(regex, substitution, text, **kwargs)

# --------------------( WASTELANDS                         )--------------------
    # The returned dictionary contains *no* entries for . Call `get_match_groups_numeric()` instead to obtain such groups.
# from sre_parse import Pattern
    # if isinstance(regex, Pattern):
    # assert isinstance(substitution, str), '"{}" not a string.'.format(
    #     substitution)
