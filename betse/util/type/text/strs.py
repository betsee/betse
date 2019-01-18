#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

#FIXME: For maintainability, shift all joiners (i.e., "join"-prefixed
#functions) into a new submodule -- say, "betse.util.type.text.strjoin".

# ....................{ IMPORTS                           }....................
import textwrap
from betse.exceptions import BetseStrException
from betse.util.type import types
from betse.util.type.types import (
    type_check, IntOrNoneTypes, IterableTypes, StrOrNoneTypes)
from textwrap import TextWrapper

# For convenience, permit callers to import the general-purpose trim() function
# from this submodule rather than the "types" submodule.
from betse.util.type.types import trim

if False: trim  # silence IDE warnings

# ....................{ SINGLETONS                        }....................
text_wrapper = TextWrapper()
'''
Singleton :class:`TextWrapper` instance with which to wrap text.

This singleton improves efficiency -- occasionally dramatically. *All* public
functions provided by the :mod:`textwrap` module implicitly instantiate
temporary :class:`TextWrapper` instances on each function call.
'''

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_if_empty(text: str, exception_message: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a human-readable
    message) if the passed string is empty.

    Parameters
    ----------
    text : str
        String to be validated.
    exception_message : optional[str]
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseStrException
        If this string is empty.
    '''

    # If this string is empty, raise an exception.
    if not text:
        # If no exception message was passed, synthesize one.
        if not exception_message:
            exception_message = 'String empty.'.format()

        # Raise this exception.
        raise BetseStrException(exception_message)


@type_check
def die_unless_prefix(
    text: str, prefix: str, exception_message: StrOrNoneTypes = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a message
    synthesized from the passed arguments) if the second passed string does
    *not* prefix the first passed string.

    Parameters
    ----------
    text : str
        String to be validated.
    prefix : str
        Prefix to test for.
    exception_message : optional[str]
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Raises
    ----------
    BetseStrException
        If this string is *not* prefixed by this prefix.
    '''

    # If this string is *NOT* prefixed by this prefix, raise an exception.
    if not is_prefix(text, prefix):
        # If no exception message was passed, synthesize one from this name.
        if not exception_message:
            exception_message = 'String "{}" not prefixed by "{}".'.format(
                text, prefix)

        # Raise this exception.
        raise BetseStrException(exception_message)

# ....................{ TESTERS ~ case                    }....................
@type_check
def is_lowercase(text: str) -> bool:
    '''
    ``True`` only if the passed string is strictly lowercase (i.e., this string
    contains at least one alphabetic character *and* all such characters are
    lowercase).

    This function ignores all non-alphabetic characters (e.g., digits,
    punctuation) in this string.
    '''

    # O.K., that's pretty sweet.
    return text.islower()


@type_check
def is_uppercase(text: str) -> bool:
    '''
    ``True`` only if the passed string is strictly uppercase (i.e., this string
    contains at least one alphabetic character *and* all such characters are
    uppercase).

    This function ignores all non-alphabetic characters (e.g., digits,
    punctuation) in this string.
    '''

    # Ditto on the sweetness.
    return text.isupper()

# ....................{ TESTERS ~ [pre|suf]fix            }....................
def is_prefix(text: str, prefix: str) -> bool:
    '''
    ``True`` only if the second passed string prefixes the first passed string.

    Parameters
    ----------
    text : str
        String to be tested.
    prefix : str
        Prefix to test for.
    '''

    return text.startswith(prefix)


@type_check
def is_suffix(text: str, suffix: str) -> bool:
    '''
    ``True`` only if the second passed string suffixes the first passed string.

    Parameters
    ----------
    text : str
        String to be tested.
    suffix : str
        Prefix to test for.
    '''

    return text.endswith(suffix)

# ....................{ GETTERS                           }....................
@type_check
def get_substr_first_index_or_none(text: str, substr: str) -> IntOrNoneTypes:
    '''
    0-based index of the passed substring in the passed string if any *or*
    ``None`` otherwise (i.e., if this string contains no such substring).

    Parameters
    ----------
    text : str
        String to be searched.
    substr : str
        Substring to search this string for.

    Returns
    ----------
    IntOrNoneTypes
        Either:

        * If this string contains this substring, the 0-based index of this
          substring in this string.
        * Else, ``None``.
    '''

    return text.index(substr) if substr in text else None

# ....................{ GETTERS ~ prefix                  }....................
@type_check
def get_prefix_preceding_char(text: str, char: str) -> str:
    '''
    **Prefix** (i.e., substring anchored at the first character) of the passed
    string preceding the first instance of the passed character in this string
    if any *or* raise an exception otherwise (i.e., if this string contains no
    such character).

    Parameters
    ----------
    text : str
        String to be searched.
    char: str
        Character to search this string for.

    Returns
    ----------
    str
        Prefix of this string preceding the first instance of this character in
        this string.

    Raises
    ----------
    ValueError
        If this string does *not* contain this character.

    See Also
    ----------
    :func:`get_prefix_preceding_char_or_none`
        Getter returning ``None`` if this string contains no such character.
    :func:`get_prefix_preceding_char_or_text`
        Getter returning ``text`` if this string contains no such character.
    '''

    # Return the prefix of this string preceding the first instance of this
    # character in this string if any *OR* raise a "ValueError" otherwise.
    return text[:text.index(char)]


@type_check
def get_prefix_preceding_char_or_none(text: str, char: str) -> StrOrNoneTypes:
    '''
    **Prefix** (i.e., substring anchored at the first character) of the passed
    string preceding the first instance of the passed character in this string
    if any *or* ``None`` otherwise (i.e., if this string contains no such
    character).

    Parameters
    ----------
    text : str
        String to be searched.
    char: str
        Character to search this string for.

    Returns
    ----------
    StrOrNoneTypes
        Either:

        * If this string contains this character, the prefix of this string
          preceding the first instance of this character in this string.
        * Else, ``None``.

    Examples
    ----------
        >>> from betse.util.type.text import strs
        >>> strs.get_prefix_preceding_char_or_none(
        ...     text='Opposition...contradiction...premonition...compromise.',
        ...     char='.')
        Opposition
        >>> strs.get_prefix_preceding_char_or_none(
        ...     text='This is an anomaly. Disabled. What is true?',
        ...     char='!')
        None
    '''

    # Return either...
    return (
        # The prefix of this string preceding the first instance of this
        # character in this string...
        text[:text.index(char)]
        # If this string contains this character *OR*...
        if char in text else
        # Nothingness.
        None
    )


@type_check
def get_prefix_preceding_char_or_text(text: str, char: str) -> str:
    '''
    **Prefix** (i.e., substring anchored at the first character) of the passed
    string preceding the first instance of the passed character in this string
    if any *or* this string as is otherwise (i.e., if this string contains no
    such character).

    Parameters
    ----------
    text : str
        String to be searched.
    char: str
        Character to search this string for.

    Returns
    ----------
    StrOrNoneTypes
        Either:

        * If this string contains this character, the prefix of this string
          preceding the first instance of this character in this string.
        * Else, this string as is.

    Examples
    ----------
        >>> from betse.util.type.text import strs
        >>> strs.get_prefix_preceding_char_or_text(
        ...     text='Seven milestones... under a watching autumn eye.',
        ...     char='.')
        Seven milestones
        >>> strs.get_prefix_preceding_char_or_text(
        ...     text='Rain's falling. Hours crawling.',
        ...     char='!')
        Rain's falling. Hours crawling.
    '''

    # Return either...
    return (
        # The prefix of this string preceding the first instance of this
        # character in this string...
        text[:text.index(char)]
        # If this string contains this character *OR*...
        if char in text else
        # This string.
        text
    )

# ....................{ ADDERS                            }....................
def add_prefix_unless_found(text: str, prefix: str) -> str:
    '''
    Prefix the passed string by the passed prefix unless the former is already
    prefixed by the latter.

    Parameters
    ----------
    text : str
        String to add this prefix to.
    prefix : str
        Prefix to add to this string.

    Returns
    ----------
    str
        The former prefixed by the latter.
    '''

    return text if is_prefix(text, prefix) else prefix + text


def add_suffix_unless_found(text: str, suffix: str) -> str:
    '''
    Suffix the passed string by the passed suffix unless the former is already
    suffixed by the latter.

    Parameters
    ----------
    text : str
        String to add this suffix to.
    suffix : str
        Suffix to add to this string.

    Returns
    ----------
    str
        The former suffixed by the latter.
    '''

    return text if is_suffix(text, suffix) else text + suffix

# ....................{ CASERS                            }....................
@type_check
def lowercase_char_first(text: str) -> str:
    '''
    Lowercase the first character of the passed string.
    '''

    return text[0].lower() + text[1:] if text else ''


@type_check
def uppercase_char_first(text: str) -> str:
    '''
    Uppercase the first character of the passed string.

    Whereas the standard :meth:`str.capitalize` method both uppercases the
    first character of this string *and* lowercases all remaining characters,
    this function *only* uppercases the first character. All remaining
    characters remain unmodified.
    '''

    return text[0].upper() + text[1:] if text else ''

# ....................{ JOINERS                           }....................
def join(*texts) -> str:
    '''
    Concatenation of the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard ``"".join((...))``
    method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='')


@type_check
def join_by_index(iterable: IterableTypes, subiterable_index: object) -> str:
    '''
    Concatenation of each string at the passed key or index of each subiterable
    of the passed iterable, with no separating delimiter.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of subiterables to be summed.
    subiterable_index : object
        Object with which to index each subiterable of this iterable. The type
        of this object *must* be a type accepted by the ``__getitem__()``
        special method of each subiterable. If each subiterable is a:

        * **Mapping** (e.g., :class:`dict`), this object *must* be hashable.
        * **Sequence** (e.g., :class:`list`, :class:`tuple`), this object
          *must* be either:

          * An integer.
          * A :func:`slice` object.

    Returns
    ----------
    str
        Concatenation of each string at this key or index of each subiterable
        of this iterable, with no separating delimiter.
    '''

    # Efficiency and simplicity combine here to form SuperHappyFunFunction.
    return ''.join(subiterable[subiterable_index] for subiterable in iterable)

# ....................{ JOINERS ~ on                      }....................
def join_on_newline(*texts) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    ``"\n".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='\n')


def join_on_dot(*texts) -> str:
    '''
    Join the passed strings with a period as the separating delimiter.

    This is a convnience function wrapping the standard
    ``".".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='.')


def join_on_space(*texts) -> str:
    '''
    Join the passed strings with a space as the separating delimiter.

    This is a convnience function wrapping the standard
    ``" ".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter=' ')


@type_check
def join_on(*texts: IterableTypes, delimiter: str) -> str:
    '''
    Join the passed strings with the passed separating delimiter.

    This is a convenience function wrapping the standard
    ``"...".join((...))`` method, whose syntax is arguably overly obfuscated.

    Parameters
    ----------
    texts : SequenceTypes
        Tuple of all strings to be joined, consisting of either:

        * One or more strings.
        * One iterable of strings.
    delimiter : str
        Substring to join each such string on.

    Returns
    ----------
    str
        String joined from the passed strings with the passed delimiter.
    '''

    # If only one object was passed and this object is a non-string iterable
    # (e.g, list, tuple), default the list of passed strings to this object.
    if len(texts) == 1 and types.is_iterable_nonstr(texts[0]):
        texts = texts[0]

    # Join these strings.
    return delimiter.join(texts)

# ....................{ JOINERS ~ as                      }....................
@type_check
def join_as(
    *texts,
    delimiter_if_two: str,
    delimiter_if_three_or_more_nonlast: str,
    delimiter_if_three_or_more_last: str
) -> str:
    '''
    Join the passed strings in a human-readable manner.

    If:

    * No strings are passed, the empty string is returned.
    * One string is passed, this string is returned as is without modification.
    * Two strings are passed, these strings are joined with the passed
      ``delimiter_if_two`` separator.
    * Three or more strings are passed:

      * All such strings except the last two are joined with the passed
        ``delimiter_if_three_or_more_nonlast`` separator.
      * The last two such strings are joined with the passed
        ``delimiter_if_three_or_more_last`` separator.

    Parameters
    ----------
    texts : tuple[str]
        Tuple of all strings to be joined.
    delimiter_if_two : str
        String separating each item of ``texts`` if ``len(texts) == 2``.
    delimiter_if_three_or_more_nonlast : str
        String separating each item *except* the last two of ``texts`` if
        ``len(texts) >= 3``.
    delimiter_if_three_or_more_last : str
        String separating the last two items of ``texts`` if
        ``len(texts) >= 3``.

    Returns
    ----------
    str
        Resulting string as described above.

    Examples
    ----------
    >>> join_as(
    ...     ('Fulgrim', 'Perturabo', 'Angron', 'Mortarion'),
    ...     delimiter_if_two=' and ',
    ...     delimiter_if_three_or_more_nonlast=', ',
    ...     delimiter_if_three_or_more_last=', and '
    ... )
    'Fulgrim, Perturabo, Angron, and Mortarion'
    '''

    # Number of passed strings.
    texts_count = len(texts)

    # If no strings are passed, return the empty string.
    if texts_count == 0:
        return ''
    # If one string is passed, return this string as is.
    elif texts_count == 1:
        return texts[0]
    # If two strings are passed, return these strings joined appropriately.
    elif texts_count == 2:
        return '{}{}{}'.format(
            texts[0], delimiter_if_two, texts[1])
    # If three or more strings are passed, return these strings joined
    # appropriately.
    else:
        # All such strings except the last two, joined appropriately.
        texts_nonlast = join_on(
            *texts[0:-2], delimiter=delimiter_if_three_or_more_nonlast)

        # The last two such strings, joined appropriately.
        texts_last = '{}{}{}'.format(
            texts[-2], delimiter_if_three_or_more_last, texts[-1])

        # Return these two substrings, joined appropriately.
        return '{}{}{}'.format(
            texts_nonlast, delimiter_if_three_or_more_nonlast, texts_last)

# ....................{ JOINERS ~ as : conjunction        }....................
@type_check
def join_as_conjunction(*texts: str) -> str:
    '''
    Conjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, and ``.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' and ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', and '
    )


@type_check
def join_as_conjunction_double_quoted(*texts: str) -> str:
    '''
    Conjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, and ``.
    '''

    # Tuple of all passed strings double-quoted. Since the "*" operator applied
    # to this tuple below requires a sequence rather than generator, this is
    # the most space-efficient available sequence (i.e., frozen tuple).
    texts_quoted = tuple(double_quote(text) for text in texts)

    # Conjunctively join these strings.
    return join_as_conjunction(*texts_quoted)

# ....................{ JOINERS ~ as : conjunction        }....................
@type_check
def join_as_disjunction(*texts: str) -> str:
    '''
    Disjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, or ``.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' or ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', or '
    )


@type_check
def join_as_disconjunction_double_quoted(*texts: str) -> str:
    '''
    Disjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, or ``.
    '''

    # Tuple of all passed strings double-quoted. See
    # join_as_conjunction_double_quoted().
    texts_quoted = tuple(double_quote(text) for text in texts)

    # Disjunctively join these strings.
    return join_as_disjunction(*texts_quoted)

# ....................{ QUOTERS                           }....................
#FIXME: We don't actually escape embedded double quotes yet. The reason why is
#that doing is subtly non-trivial; for safety, we don't want to re-escape
#already escaped double quotes in this string. We only want to escape
#currently unescaped double quotes in this string. To do so, we'll need to
#leverage regex-based substitution with negative lookbehind preventing existing
#substrings matching '\\"' from being subject to substitution.

def double_quote(text: str) -> str:
    '''
    Double-quote the passed string in a human-readable manner.

    Specifically (in order):

    #. All prefixing and suffixing whitespace is removed from this string.
    #. If this string is either double- or single-quoted, these quotes are
       removed.
    #. All other double quotes in this string are escaped (e.g., replacing each
       ``"`` character with ``\\"``).
    #. The resulting string is double-quoted and returned.
    '''

    # Remove all prefixing and suffixing whitespace from this string.
    text = remove_whitespace_presuffix(text)

    # If this string is already either double- or single-quoted, strip these
    # delimiting quotes for simplicity.
    if ((text[0] == '"' and text[-1] == '"') or
        (text[0] == "'" and text[-1] == "'")):
        text = text[1:-1]

    # Double quote this string.
    return '"{}"'.format(text)

# ....................{ REMOVERS ~ prefix                 }....................
@type_check
def remove_prefix(
    text: str, prefix: str, exception_message: str = None) -> str:
    '''
    Passed string with the passed prefix removed if present *or* raise an
    exception with the passed message otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    prefix : str
        Prefix to remove from this string.
    exception_message : optional[str]
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Returns
    ----------
    str
        This string truncated as detailed above.
    '''

    # If this string is *NOT* prefixed by this prefix, raise this exception.
    die_unless_prefix(text, prefix, exception_message)

    # Else, return a new string with this prefix removed from this string.
    return remove_prefix_if_found(text, prefix)


@type_check
def remove_prefix_if_found(text: str, prefix: str) -> str:
    '''
    Passed string with the passed prefix removed if present *or* the passed
    string as is otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    prefix : str
        Prefix to remove from this string.

    Returns
    ----------
    str
        This string truncated as detailed above.
    '''

    return text[len(prefix):] if is_prefix(text, prefix) else text

# ....................{ REMOVERS ~ suffix                 }....................
@type_check
def remove_suffix_if_found(text: str, suffix: str) -> str:
    '''
    Passed string with the passed suffix removed if present *or* the passed
    string as is otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    suffix : str
        Suffix to remove from this string.

    Returns
    ----------
    str
        This string truncated as detailed above.
    '''

    # There exists a special case *NOT* present in remove_prefix(). If this
    # suffix is empty, "string[:-0]" is also incorrectly empty. Avoid returning
    # the empty string in this case by explicitly testing for emptiness.
    return text[:-len(suffix)] if suffix and is_suffix(text, suffix) else text


@type_check
def remove_suffix_prefixed(text: str, suffix_prefix: str) -> str:
    '''
    Passed string with *all* characters including and following the first
    instance of the passed substring if present removed *or* the passed string
    as is otherwise.

    Specifically, this functions returns:

    * If this substring is the empty string, the empty string.
    * Else if this string contains no such substrings, this string as is.
    * Else, only the prefix of this string preceding the first such substring
      in this string.

    Parameters
    ----------
    text : str
        String to be inspected.
    suffix_prefix : str
        Substring of this string to begin removing characters at.

    Returns
    ----------
    str
        This string truncated as detailed above.

    See Also
    ----------
    :func:`remove_suffix_prefixed`
        Function replacing rather than removing this suffix.
    '''

    # If:
    #
    # * This suffix prefix is non-empty, return the substring prefixing this
    #   suffix prefix in this string if any *OR* this string as is otherwise.
    #   Fortuitously, this is exactly the first value returned by the low-level
    #   str.partition() method underlying this implementation.
    # * This suffix prefix is the empty string, return the empty string. This
    #   edge case must be explicitly tested, as str.partition() raises this
    #   exception in this case: "ValueError: empty separator".
    return text.partition(suffix_prefix)[0] if suffix_prefix else ''

# ....................{ REMOVERS ~ space                  }....................
@type_check
def remove_whitespace(text: str) -> str:
    '''
    Passed string with *all* whitespace removed, including all prefixing and
    suffixing whitespace and whitespace interspersed through this string.

    Parameters
    ----------
    text : str
        String to remove all whitespace from.

    Returns
    ----------
    str
        String with all whitespace removed.

    See Also
    ----------
    https://stackoverflow.com/a/8270124/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # While a regular expression-based implementation is also feasible, the
    # current approach of splitting this string on whitespace into substrings
    # containing no whitespace and then concatenating these substrings is
    # substantially more Pythonic (and presumably efficient).
    return ''.join(text.split())


@type_check
def remove_whitespace_presuffix(text: str) -> str:
    '''
    Passed string with all prefixing and suffixing whitespace removed.
    '''

    return text.strip()


@type_check
def remove_whitespace_suffix(text: str) -> str:
    '''
    Passed string with all suffixing (but *not* prefixing) whitespace removed.
    '''

    return text.rstrip()

# ....................{ REMOVERS ~ space : newline        }....................
@type_check
def remove_newlines_suffix(text: str) -> str:
    '''
    Passed string with all suffixing (but *not* prefixing) newlines removed.
    '''

    return text.rstrip('\n')

# ....................{ REPLACERS                         }....................
@type_check
def replace_substrs(text: str, substr: str, replacement: str) -> str:
    '''
    Passed string with all instances of the passed substring replaced by the
    passed replacement substring if any *or* the passed string as is otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    substr : str
        Substring to replace all instances of in this string.
    replacement : str
        Substring to replace all instances of the passed substring.

    Returns
    ----------
    str
        This string with all instances of this substring replaced by this
        replacement substring.

    See Also
    ----------
    :func:`betse.util.type.text.regexes.replace_substrs`
        Equivalent regular expression-based replacer.
    '''

    # ...thas just how we roll.
    return text.replace(substr, replacement)


@type_check
def replace_substr_first(text: str, substr: str, replacement: str) -> str:
    '''
    Passed string with the first instance of the passed substring replaced by
    the passed replacement substring if any *or* the passed string as is
    otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    substr : str
        Substring to replace the first instance of in this string.
    replacement : str
        Substring to replace the first instance of the passed substring.

    Returns
    ----------
    str
        This string with the first instance of this substring replaced by this
        replacement substring.
    '''

    # ...thas just how we roll.
    return text.replace(substr, replacement, 1)

# ....................{ REPLACERS ~ suffix                }....................
@type_check
def replace_suffix_prefixed(
    text: str, suffix_prefix: str, replacement: str) -> str:
    '''
    Passed string with *all* characters including and following the first
    instance of the passed substring if present replaced by the passed
    replacement substring *or* this string returned as is otherwise.

    Parameters
    ----------
    text : str
        String to be inspected.
    suffix_prefix : str
        Substring of this string to begin replacing characters at.
    replacement : str
        Substring to replace the suffix of this string prefixed by this suffix
        prefix by.

    Returns
    ----------
    str
        This string with all characters including and following the first
        instance of the passed substring replaced by this replacement
        substring.

    See Also
    ----------
    :func:`remove_suffix_prefixed`
        Function removing rather than replacing this suffix.
    '''

    # 0-based index of the first instance of this suffix prefix in this string
    # if present *OR* "None" otherwise.
    suffix_prefix_first_index = get_substr_first_index_or_none(
        text=text, substr=suffix_prefix)

    # Return...
    return (
        # If this string contains this suffix prefix, the prefix of this string
        # preceding this suffix prefix appended by this replacement substring.
        text[:suffix_prefix_first_index] + replacement
        if suffix_prefix_first_index is not None else
        # Else, this string as is.
        text
    )

# ....................{ TRUNCATERS                        }....................
@type_check
def truncate(
    text: str,
    max_len: int = 80,
    suffix_prefix: StrOrNoneTypes = None,
    replacement: str = '...',
) -> str:
    '''
    Passed string truncated to the passed suffix prefix and maximum length if
    applicable *or* this string as is otherwise.

    Specifically, this function applies the following operations (in order):

    #. If a suffix prefix is passed, the suffix of the passed string prefixed
       by this prefix is replaced by the passed replacement.
    #. If the resulting string still exceeds the passed maximum length, the
       suffix of this string exceeding this length is again replaced by the
       passed replacement.
    #. The resulting string is returned.

    Parameters
    ----------
    text : str
        String to be truncated.
    max_len : int
        Maximum number of characters to truncate this string to. Defaults to
        the standard UNIX terminal line length (i.e., 80).
    suffix_prefix : StrOrNoneTypes
        Substring of this string to begin prematurely truncating characters at
        if this string contains this substring, regardless of whether this
        string exceeds this maximum length. This substring imposes a "barrier."
        Defaults to ``None``, in which case no such barrier is imposed.
    replacement : str
        Substring to replace the truncated portion of this string with.
        Defaults to an ASCII ellipses (i.e., ``...``).

    Returns
    ----------
    str
        This string truncated to this maximum length, as detailed above.
    '''

    # If passed a suffix prefix, replace the suffix of this string prefixed
    # by this prefix by this replacement.
    if suffix_prefix is not None:
        text = replace_suffix_prefixed(
            text=text, suffix_prefix=suffix_prefix, replacement=replacement)

    # If this string does *NOT* exceed this maximum, return this string as is.
    if len(text) <= max_len:
        return text
    # Else, this string exceeds this maximum. In this case...
    else:
        # Number of characters to truncate from the end of this string.
        # Dismantled, this is:
        #
        # * "len(text) - max_len", the number of characters that this string
        #   exceeds this maximum length by.
        # * "... + len(replacement)", truncating an additional number of
        #   characters equal to the length of this replacement so as to make
        #   sufficient room for this replacement at the end of this string
        #   without exceeding this maximum length.
        #
        # Note that this number is guaranteed to be non-negative (i.e., greater
        # than zero), as "len(text) > max_len" and "len(replacement) >= 0".
        truncate_chars_count = len(text) - max_len + len(replacement)

        # If more characters are to be truncated from this string than exist in
        # this string, then it can be shown by trivial algebraic equivalency
        # that "len(replacement) > max_len" (i.e., the length of the
        # replacement substring exceeds the maximum length). In this uncommon
        # edge case, return the replacement truncated to this maximum.
        if truncate_chars_count > len(text):
            return replacement[:max_len]
        # Else, fewer characters are to be truncated from this string than
        # exist in this string. This is the common case.

        # Return this string truncated to this number of characters appended by
        # this replacement substring.
        return text[:-truncate_chars_count] + replacement

# ....................{ WRAPPERS                          }....................
@type_check
def wrap_lines(lines: IterableTypes, **kwargs) -> str:
    '''
    Wrap the passed iterable of lines to the passed line width, prefixing each
    resulting wrapped line by the passed line prefix.

    See Also
    ----------
    :func:`wrap`
        Further details.
    '''

    return wrap(join(lines), **kwargs)


@type_check
def wrap(
    text: str,
    text_wrapper = textwrap,
    line_prefix: str = '',
    **kwargs
) -> str:
    '''
    Wrap the passed text to the passed line width, prefixing each resulting
    wrapped line by the passed line prefix.

    This function accepts the following keyword arguments in addition to those
    accepted by the :func:`textwrap.wrap` function:

    * ``text_wrapper`, the object on which to call the :func:`wrap` function or
      method. For safety, this defaults to the :mod:`textwrap` module.
    * ``line_prefix``, the substring prefixing each wrapped output line.

    See Also
    ----------
    https://docs.python.org/3/library/textwrap.html
        For further details on keyword arguments.
    '''

    assert hasattr(text_wrapper, 'wrap'), (
        'Object "{}" not a text wrapper '
        '(i.e., has no wrap() callable).'.format(text_wrapper))

    # wrap() function or method to be called.
    wrap_callable = getattr(text_wrapper, 'wrap')
    assert callable(wrap_callable), (
        'Text wrapper "{}" attribute "wrap" not callable.'.format(
            text_wrapper))

    # If passed a nonempty line prefix, add appropriate keyword arguments.
    if line_prefix:
        kwargs['initial_indent'] = line_prefix
        kwargs['subsequent_indent'] = line_prefix

    # For reliability, call the wrap() function of module "textwrap" rather than
    # the singleton object "text_wrapper". The former internally instantiates a
    # new instance of class "TextWrapper" and hence resets all wrapping
    # attributes to sensible defaults, whereas the latter reuses existing such
    # attributes -- which may no longer retain sensible defaults.
    return join_on_newline(wrap_callable(text, **kwargs))

# ....................{ WRAPPERS ~ un                     }....................
@type_check
def unwrap(text: str) -> str:
    '''
    Convert the passed possibly **multiline string** (i.e., string containing
    one or more newlines) to a **single-line string** (i.e., string containing
    no newlines).
    '''

    # Nice one, stdlib.
    return text.replace('\n', ' ')

# ....................{ (IN|DE)DENTERS                    }....................
def dedent(*texts) -> str:
    '''
    Remove all indentation shared in common by all lines of all passed strings.
    '''

    return textwrap.dedent(*texts)
