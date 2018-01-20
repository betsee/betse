#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **regex** (i.e., Python-compatible regular expression) facilities.
'''

# ....................{ IMPORTS                            }....................
import re
from betse.exceptions import BetseRegexException
from betse.util.type.types import (
    type_check,
    CallableOrStrTypes,
    IterableTypes,
    MappingType,
    RegexCompiledType,
    RegexMatchType,
    RegexMatchOrNoneTypes,
    RegexTypes,
    SequenceTypes,
    SequenceOrNoneTypes,
)

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

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_replace_substrs_line(
    text: str,
    regex: RegexTypes,
    replacement: CallableOrStrTypes,
    **kwargs
) -> str:
    '''
    Passed subject string with all substrings matching the passed regular
    expression replaced by the passed substitution in a line-oriented manner if
    this subject string contains at least one such substring *or* raise an
    exception otherwise.

    Raises
    ----------
    BetseRegexException
        If this subject string fails to match this regular expression.

    See Also
    ----------
    :func:`replace_substrs_line`
        Further details.
    '''

    # If this subject string fails to match this regex, raise an exception.
    #
    # Sadly, the re.sub() function called by the replace_substrs_line() function
    # called below fails to return metadata detailing the number of replacements
    # performed or even whether any replacements were performed. Hence, this
    # validation *MUST* be performed manually beforehand.
    if not is_match_line(text=text, regex=regex, **kwargs):
        raise BetseRegexException(
            'Subject string "{}" not matched by '
            'regular expression "{}".'.format(text, regex))

    # Else, this subject string matches this regex. In this case, perform this
    # search-and-replacement and return the result.
    return replace_substrs_line(
        text=text,
        regex=regex,
        replacement=replacement,
        **kwargs
    )

# ....................{ TESTERS                            }....................
def is_match(text: str, regex: RegexTypes, **kwargs) -> bool:
    '''
    ``True`` only if zero or more characters anchored to the beginning of the
    passed string match the passed regular expression.

    Parameters
    ----------
    text : str
        String to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    bool
        ``True`` only if this string matches this regular expression.
    '''

    return get_match_or_none(text, regex, **kwargs) is not None


def is_match_line(text: str, regex: RegexTypes, **kwargs) -> bool:
    '''
    ``True`` only if one or more lines of the passed subject string match the
    passed regular expression.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.search` function.

    Returns
    ----------
    bool
        ``True`` only if this string matches this regular expression.

    See Also
    ----------
    :func:`get_match_line_first_or_none`
        Further details.
    '''

    return get_match_line_first_or_none(text, regex, **kwargs) is not None

# ....................{ MATCHERS ~ group : named           }....................
def get_match_groups_named(
    text: str, regex: RegexTypes, **kwargs) -> MappingType:
    '''
    Dictionary mapping explicitly named groups to substrings matched anchored to
    the beginning of the passed string against the passed regular expression if
    a match exists or raise an exception otherwise.

    For each capture group of the form ``(?P<group_name>...)`` in this regular
    expression, this dictionary contains a key-value pair whose:

    * Key is this named group's ``group_name``.
    * Value is either:
      * ``None`` if this group is unmatched by this string.
      * The matched substring otherwise.

    Unnamed (i.e., only numbered) groups are ignored and hence excluded from
    this dictionary, regardless of whether any of these groups matched. If
    undesirable, call the :func:`get_match_groups_numbered` function instead.

    Caveats
    ----------
    Python does *not* support the increasingly standardized named group form
    ``(?<group_name>...)`` -- only the Python-specific form
    ``(?P<group_name>...)``.

    Parameters
    ----------
    text : str
        String to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    MappingType
        Dictionary mapping matched named groups.

    Raises
    ----------
    BetseRegexException
        If this string does *not* match this expression.

    See Also
    ----------
    :func:`get_match_or_none`
        Further details on regular expressions and keyword arguments.
    '''

    return get_match(text, regex, **kwargs).groupdict()

# ....................{ MATCHERS ~ group : numbered        }....................
def get_match_groups_numbered(
    text: str, regex: RegexTypes, **kwargs) -> SequenceTypes:
    '''
    Sequence of all groups matched anchored to the beginning of the passed
    subject string against the passed regular expression (ordered by the
    left-to-right lexical position at which each such group was matched) if any
    *or* raise an exception otherwise.

    Unmatched groups will have the value `None`.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    SequenceTypes
        Sequence of all matched groups.

    Raises
    ----------
    BetseRegexException
        If this string does *not* match this expression.

    See Also
    ----------
    :func:`get_match_or_none`
        Further details.
    '''

    return get_match(text, regex, **kwargs).groups()


def get_match_groups_numbered_or_none(
    text: str, regex: RegexTypes, **kwargs) -> SequenceOrNoneTypes:
    '''
    Sequence of all groups matched anchored to the beginning of the passed
    subject string against the passed regular expression (ordered by the
    left-to-right lexical position at which each such group was matched) if any
    *or* ``None`` otherwise.

    Unmatched groups will have the value ``None``.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    SequenceOrNoneTypes
        Either:
        * If this string matches this regular expression, the sequence of all
          groups matched from this string.
        * Else, ``None``.

    See Also
    ----------
    :func:`get_match_or_none`
        Further details.
    '''

    match = get_match_or_none(text, regex, **kwargs)
    return match.groups() if match is not None else None

# ....................{ MATCHERS ~ full : first            }....................
def get_match_full_first_or_none(
    text: str, regex: RegexTypes, **kwargs) -> SequenceOrNoneTypes:
    '''
    First complete substring matched anchored to the beginning of the passed
    subject string against the passed regular expression if any *or* ``None``
    otherwise.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    SequenceOrNoneTypes
        Either:
        * If this string matches this regular expression, the first complete
          substring matched from this string.
        * Else, ``None``.

    See Also
    ----------
    :func:`get_match_or_none`
        Further details.
    '''

    match = get_match_or_none(text, regex, **kwargs)
    return match.group(0) if match is not None else None

# ....................{ MATCHERS ~ obj                     }....................
def get_match(text: str, regex: RegexTypes, **kwargs) -> RegexMatchType:
    '''
    Match object obtained by matching zero or more characters anchored to the
    beginning of the passed subject string against the passed regular expression
    if any match exists *or* raise an exception otherwise.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.match` function.

    Returns
    ----------
    RegexMatchType
        This match object.

    Raises
    ----------
    BetseRegexException
        If this string does *not* match this expression.

    See Also
    ----------
    :func:`get_match_or_none`
        Further details on calling conventions.
    '''

    # Object matching this string against this expression.
    match = get_match_or_none(text, regex, **kwargs)

    # If no match was found, convert the non-fatal "None" returned by the
    # re.match() function into a fatal exception. By design, no callables in the
    # "re" module raise exceptions.
    if match is None:
        raise BetseRegexException(
            'Subject string "{}" not matched by '
            'regular expression "{}".'.format(text, regex))

    # Return this match.
    return match


@type_check
def get_match_or_none(
    text: str, regex: RegexTypes, **kwargs) -> RegexMatchOrNoneTypes:
    '''
    Match object obtained by matching zero or more characters anchored to the
    beginning of the passed subject string against the passed regular expression
    if any match exists *or* ``None`` otherwise.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as
    :func:`re.match`.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.DOTALL`, forcing the ``.`` special character to match any
      character including newline. By default, this character matches any
      character excluding newline. The former is almost always preferable.

    Returns
    ----------
    RegexMatchOrNoneTypes
        This match object if a match exists *or* ``None`` otherwise.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags.
    _init_kwargs_flags_nonline(regex, kwargs)

    # Return the only object matching this string against this expression.
    return re.match(regex, text, **kwargs)

# ....................{ MATCHERS ~ obj : first             }....................
@type_check
def get_match_first_or_none(
    text: str, regex: RegexTypes, **kwargs) -> RegexMatchOrNoneTypes:
    '''
    First match object obtained by matching the passed subject string against
    the passed regular expression if any match exists *or* ``None`` otherwise.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as
    :func:`re.match`.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.DOTALL`, forcing the ``.`` special character to match any
      character including newline. By default, this character matches any
      character excluding newline. The former is almost always preferable.

    Returns
    ----------
    RegexMatchOrNoneTypes
        This match object if a match exists *or* ``None`` otherwise.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.match
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags.
    _init_kwargs_flags_nonline(regex, kwargs)

    # Return the first object matching this string against this expression.
    return re.search(regex, text, **kwargs)

# ....................{ MATCHERS ~ obj : line : first      }....................
@type_check
def get_match_line_first_or_none(
    text: str, regex: RegexTypes, **kwargs) -> RegexMatchOrNoneTypes:
    '''
    First match object obtained by matching the passed subject string against
    the passed regular expression in a line-oriented manner if any such match
    exists *or* ``None`` otherwise.

    To ensure that only single lines are matched, this regular expression should
    typically be prefixed by the ``^`` special character and/or suffixed by the
    ``$`` special character, anchoring matches to the start and/or end of this
    subject string as a whole and each line of this string.

    Parameters
    ----------
    text : str
        Subject string to match.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.search` function.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.MULTILINE`, forcing the ``^`` and ``$`` special characters to
      match both at the start and end of this string *and* at the start and end
      of each line of this string. By default, these characters match only at
      the former. Line-oriented matching requires both, however.

    Returns
    ----------
    RegexMatchOrNoneTypes
        This match object if a match exists *or* ``None`` otherwise.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.search
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags for line-oriented matching.
    _init_kwargs_flags_line(regex, kwargs)

    # Return the first object matching this string against this expression.
    return re.search(regex, text, **kwargs)

# ....................{ ITERATORS                          }....................
@type_check
def iter_matches(text: str, regex: RegexTypes, **kwargs) -> IterableTypes:
    '''
    Generator iteratively yielding each non-overlapping match at any position of
    the passed string against the passed regular expression as a match object.

    If no such match exists, this function successfully returns the empty
    generator rather than raising a fatal exception.

    Parameters
    ----------
    text : str
        Subject string to match on.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.finditer` function.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.DOTALL`, forcing the ``.`` special character to match any
      character including newline. By default, this character matches any
      character excluding newline. The former is almost always preferable.

    Returns
    ----------
    IterableTypes
        Generator yielding match objects (i.e., instances of `re.SRE_Match`).

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.search
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags.
    _init_kwargs_flags_nonline(regex, kwargs)

    # Return this generator.
    return re.finditer(regex, text, **kwargs)


@type_check
def iter_matches_line(text: str, regex: RegexTypes, **kwargs) -> IterableTypes:
    '''
    Generator iteratively yielding each non-overlapping match at any position of
    the passed string against the passed regular expression in a line-oriented
    manner as a match object.

    If no such match exists, this function successfully returns the empty
    generator rather than raising a fatal exception.

    To ensure that only single lines are matched, this regular expression should
    typically be prefixed by the `^` special character and/or suffixed by the
    `$` special character, thus anchoring matches to the start and/or end of
    this string as a whole and each line of this string.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.MULTILINE`, forcing the ``^`` and ``$`` special characters to
      match both at the start and end of this string *and* at the start and end
      of each line of this string. By default, these characters match only at
      the former. Line-oriented matching requires both, however.

    Parameters
    ----------
    text : str
        Subject string to match on.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as the
    :func:`re.finditer` function.

    Returns
    ----------
    IterableTypes
        Generator yielding match objects (i.e., instances of `re.SRE_Match`).
    '''

    # Sanitize the passed match flags for line-oriented matching.
    _init_kwargs_flags_line(regex, kwargs)

    # Return this generator.
    return re.finditer(regex, text, **kwargs)

# ....................{ REMOVERS                           }....................
def remove_substrs(text: str, regex: RegexTypes, **kwargs) -> str:
    '''
    Remove all substrings in the passed string matching the passed regular
    expression.

    Parameters
    ----------
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.

    This function accepts the same optional keyword arguments as
    :func:`re.sub`.

    Returns
    ----------
    str
        Passed string containing no substrings matching such regular expression.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''

    return replace_substrs(text=text, regex=regex, replacement='', **kwargs)

# ....................{ REPLACERS                          }....................
@type_check
def replace_substrs(
    text: str,
    regex: RegexTypes,
    replacement: CallableOrStrTypes,
    **kwargs
) -> str:
    '''
    Passed subject string with all substrings matching the passed regular
    expression replaced by the passed substitution in a non-line-oriented
    manner.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.DOTALL`, forcing the ``.`` special character to match any
      character including newline. By default, this character matches any
      character excluding newline. The former is almost always preferable.

    Parameters
    ----------
    text : str
        Subject string to replace these substrings of.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.
    replacement : CallableOrStrTypes
        Substitution to be performed, either a:
        * String.
        * Callable (e.g., function, lambda, method).

    This function accepts the same optional keyword arguments as the
    :func:`re.sub` function.

    Returns
    ----------
    str
        Passed string with all substrings matching this regular expression
        globally replaced with this substitution.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags.
    _init_kwargs_flags_nonline(regex, kwargs)

    # Substitute, if you please.
    return re.sub(regex, replacement, text, **kwargs)

# ....................{ REPLACERS : line                   }....................
@type_check
def replace_substrs_line(
    text: str,
    regex: RegexTypes,
    replacement: CallableOrStrTypes,
    **kwargs
) -> str:
    '''
    Passed subject string with all substrings matching the passed regular
    expression replaced by the passed substitution in a line-oriented manner.

    If this subject string contains no such substrings, this function silently
    reduces to a noop. If undesirable, consider calling the stricter
    :func:`die_unless_replace_substrs_line` function instead.

    Match Flags
    ----------
    For convenience, the following match flags will be enabled by default:

    * :data:`re.MULTILINE`, forcing the ``^`` and ``$`` special characters to
      match both at the start and end of this string *and* at the start and end
      of each line of this string. By default, these characters match only at
      the former. Line-oriented matching requires both, however.

    Parameters
    ----------
    text : str
        Subject string to replace these substrings of.
    regex : RegexTypes
        Regular expression to be matched. This object should be either of type:
        * :class:`str`, signifying an uncompiled regular expression.
        * :class:`Pattern`, signifying a compiled regular expression object.
    replacement : CallableOrStrTypes
        Substitution to be performed, either a:
        * String.
        * Callable (e.g., function, lambda, method).

    This function accepts the same optional keyword arguments as the
    :func:`re.sub` function.

    Returns
    ----------
    str
        Passed string with all substrings matching this regular expression
        globally replaced with this substitution.

    See Also
    ----------
    https://docs.python.org/3/library/re.html#re.sub
        Further details on regular expressions and keyword arguments.
    '''

    # Sanitize the passed match flags in a line-oriented manner.
    _init_kwargs_flags_line(regex, kwargs)

    # Substitute, if you please.
    return re.sub(regex, replacement, text, **kwargs)

# ....................{ COMPILERS                          }....................
@type_check
def compile_regex(regex: str, **kwargs) -> RegexCompiledType:
    '''
    Compile the passed uncompiled regular expression.

    Parameters
    ----------
    regex : str
        Uncompiled regular expression to be compiled.

    All remaining keyword parameters are passed as is to the :func:`re.compile`
    function.

    Returns
    ----------
    RegexCompiledType
        This compiled regular expression.
    '''

    # Return this regular expression compiled.
    return re.compile(regex, **kwargs)

# ....................{ SUBSTITUTERS                       }....................
@type_check
def _init_kwargs_flags_nonline(regex: RegexTypes, kwargs: MappingType) -> None:
    '''
    Sanitize the list of match flags in the passed dictionary for
    non-line-oriented matching.

    Specifically, this function adds the :data:`re.DOTALL` flag to the integer
    value of the ``flags`` key of this dictionary (defaulting to zero if
    currently unset).
    '''

    # If this regular expression is already compiled, reduce to a noop. Why?
    # Because flags *CANNOT* be respecified after the compilation phase.
    if isinstance(regex, RegexCompiledType):
        return

    # Else, this regular expression is uncompiled. In this case, these flags are
    # safely modifiable as required.
    kwargs['flags'] = kwargs.get('flags', 0) | re.DOTALL


@type_check
def _init_kwargs_flags_line(regex: RegexTypes, kwargs: MappingType) -> None:
    '''
    Sanitize the list of match flags in the passed dictionary for line-oriented
    matching.

    Specifically, this function adds the :data:`re.MULTILINE` flag to the
    integer value of the ``flags`` key of this dictionary (defaulting to zero if
    currently unset).
    '''

    # If this regular expression is already compiled, reduce to a noop. Why?
    # Because flags *CANNOT* be respecified after the compilation phase.
    if isinstance(regex, RegexCompiledType):
        return

    # Else, this regular expression is uncompiled. In this case, these flags are
    # safely modifiable as required.
    kwargs['flags'] = kwargs.get('flags', 0) | re.MULTILINE
