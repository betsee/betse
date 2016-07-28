#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **regex** (i.e., Python-compatible regular expression) facilities.
'''

# ....................{ IMPORTS                            }....................
import re
from betse.exceptions import BetseRegexException
from betse.util.type.types import type_check, CallableTypes
from sre_parse import Pattern

# ....................{ FLAGS                              }....................
FLAG_MULTILINE = re.MULTILINE
'''
When specified, the pattern character:

* `^` matches both at the beginning of the subject string _and_ at the beginning
  of each line (immediately following each newline) of this string.
* `$` matches both at the end of the subject string _and_ at the end of each
  line (immediately preceding each newline) of this string.

By default:

* `^` matches only at the beginning of the subject string.
* `$` matches only at the end of the subject string and immediately before the
  newline (if any) at the end of this string.
'''

# ....................{ TESTERS                            }....................
def is_match(text: str, regex: (str, Pattern), **kwargs) -> bool:
    '''
    `True` only if the passed string matches the passed regular expression.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Returns
    ----------
    bool
        `True` only if this string matches this regular expression.
    '''

    return get_match_if_any(text, regex, **kwargs) is not None

# ....................{ MATCHERS ~ group : named           }....................
def get_match_groups_named(text: str, regex: (str, Pattern), **kwargs) -> list:
    '''
    Get the dictionary mapping explicitly named groups to substrings matched
    from the passed string against the passed regular expression if a match
    exists or raise an exception otherwise.

    Unmatched groups will have the value `None`. Unnamed (i.e., only numbered)
    groups will be ignored and hence absent from this dictionary, regardless of
    whether any such group matched. If this is undesirable, consider calling
    `get_match_groups_numbered()` instead.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Returns
    ----------
    dict
        Dictionary mapping matched named groups.

    Raises
    ----------
    BetseRegexException
        If this string does _not_ match this expression.

    See Also
    ----------
    get_match_if_any
        Further details on regular expressions and keyword arguments.
    '''

    return get_match(text, regex, **kwargs).groupdict()

# ....................{ MATCHERS ~ group : numbered        }....................
def get_match_groups_numbered(
    text: str, regex: (str, Pattern), **kwargs) -> list:
    '''
    Get the list of all groups matched from the passed string against the passed
    regular expression, ordered by the left-to-right lexical position at which
    each such group was matched, if a match exists or raise an exception
    otherwise.

    Unmatched groups will have the value `None`.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Returns
    ----------
    list
        List of matched groups.

    Raises
    ----------
    BetseRegexException
        If this string does _not_ match this expression.

    See Also
    ----------
    get_match_if_any
        Further details on regular expressions and keyword arguments.
    '''

    return get_match(text, regex, **kwargs).groups()


def get_match_groups_numbered_if_any(
    text: str, regex: (str, Pattern), **kwargs) -> list:
    '''
    Get the list of all groups matched from the passed string against the passed
    regular expression, ordered by the left-to-right lexical position at which
    each such group was matched, if a match exists or `None` otherwise.

    Unmatched groups will have the value `None`.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Returns
    ----------
    list
        List of matched groups if any or `None` otherwise.

    See Also
    ----------
    get_match_if_any
        Further details on regular expressions and keyword arguments.
    '''

    match = get_match_if_any(text, regex, **kwargs)
    return match.groups() if match is not None else None

# ....................{ MATCHERS ~ object                  }....................
def get_match(text: str, regex: (str, Pattern), **kwargs) -> 'SRE_Match':
    '''
    Get the match object obtained by matching the passed string against the
    passed regular expression if any match exists _or_ raise an exception
    otherwise.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Returns
    ----------
    SRE_Match
        The match object.

    Raises
    ----------
    BetseRegexException
        If this string does _not_ match this expression.

    See Also
    ----------
    get_match_if_any
        Further details on calling conventions.
    '''

    # Match group of this string against this expression.
    match = get_match_if_any(text, regex, **kwargs)

    # If no match was found, convert the non-fatal "None" returned by re.match()
    # into a fatal exception. By design, no callables of the standard re module
    # raise exceptions.
    if match is None:
        raise BetseRegexException(
            'Subject string "{}" not matched by '
            'regular expression "{}".'.format(text, regex))

    return match


@type_check
def get_match_if_any(text: str, regex: (str, Pattern), **kwargs):
    '''
    Get the match object obtained by matching the passed string against the
    passed regular expression if any match exists _or_ `None` otherwise.

    Parameters
    ----------
    text : str
        String to match.
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.match()`.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * `re.DOTALL`, forcing the `.` special character to match any character
      including newline. By default, this character matches any character
      excluding newline. The former is almost always preferable, however.

    Returns
    ----------
    re.SRE_Match, None
        The match object if a match exists _or_ `None` otherwise.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''

    # Enable the following match flags by default.
    kwargs['flags'] = kwargs.get('flags', 0) | re.DOTALL

    # Match group of this string against this expression.
    return re.match(regex, text, **kwargs)

# ....................{ REPLACERS                          }....................
def remove_substrings(text: str, regex: (str, Pattern), **kwargs) -> str:
    '''
    Remove all substrings in the passed string matching the passed regular
    expression.

    Parameters
    ----------
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as `re.sub()`.

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
@type_check
def substitute_substrings(
    text: str,
    regex: (str, Pattern),
    substitution: (str,) + CallableTypes,
    **kwargs
) -> str:
    '''
    Replace all substrings in the passed string matching the passed regular
    expression if any with the passed substitution _or_ noop otherwise.

    Parameters
    ----------
    regex : str, Pattern
        Regular expression to be matched. This object should be either of type:
        * `str`, signifying an uncompiled regular expression.
        * `Pattern`, signifying a compiled regular expression object.
    substitution : str, func
        Substitution to be performed. This should be either a string or
        callable (e.g., function, method).

    This function accepts the same optional keyword arguments as `re.sub()`.

    Returns
    ----------
    str
        Passed string with all substrings matching this regular expression
        replaced with this substitution.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''

    # Enable the following substitution flags by default.
    kwargs['flags'] = kwargs.get('flags', 0) | re.DOTALL

    # Substitute, if you please.
    return re.sub(regex, substitution, text, **kwargs)
