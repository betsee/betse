#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

# ....................{ IMPORTS                            }....................
import textwrap
from betse.exceptions import BetseStrException
from betse.util.type import types
from betse.util.type.types import type_check, IterableTypes, StrOrNoneTypes
from textwrap import TextWrapper

# For convenience, permit callers to import the general-purpose trim() function
# from this submodule rather than the "types" submodule.
from betse.util.type.types import trim
if False: trim  # silence IDE warnings

# ....................{ SINGLETONS                         }....................
text_wrapper = TextWrapper()
'''
Singleton :class:`TextWrapper` instance with which to wrap text.

Such singleton improves efficiency -- occasionally dramatically. *All* public
functions provided by the :mod:`textwrap` module implicitly instantiate
temporary :class:`TextWrapper` instances on each call to such functions.
'''

# ....................{ EXCEPTIONS                         }....................
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
            exception_message = 'Text empty.'.format()

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
            exception_message = 'Text "{}" not prefixed by "{}".'.format(
                text, prefix)

        # Raise this exception.
        raise BetseStrException(exception_message)

# ....................{ TESTERS ~ case                     }....................
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

# ....................{ TESTERS ~ [pre|suf]fix             }....................
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

# ....................{ ADDERS                             }....................
def add_prefix_unless_found(text: str, prefix: str) -> str:
    '''
    Prefix the passed string by the passed prefix unless such string is already
    prefixed by such prefix.
    '''

    return text if is_prefix(text, prefix) else prefix + text


def add_suffix_unless_found(text: str, suffix: str) -> str:
    '''
    Suffix the passed string by the passed suffix unless such string is already
    suffixed by such suffix.
    '''

    return text if is_suffix(text, suffix) else text + suffix

# ....................{ JOINERS ~ on                       }....................
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
        Concatenation of each string at this key or index of each subiterable of
        this iterable, with no separating delimiter.
    '''

    # Efficiency and simplicity combine here to form SuperHappyFunFunction.
    return ''.join(subiterable[subiterable_index] for subiterable in iterable)

# ....................{ JOINERS ~ on                       }....................
def join_on_newline(*texts) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    ``"\n".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='\n')


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

# ....................{ JOINERS ~ as                       }....................
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
      `delimiter_if_two` separator.
    * Three or more strings are passed:
      * All such strings except the last two are joined with the passed
        `delimiter_if_three_or_more_nonlast` separator.
      * The last two such strings are joined with the passed
        `delimiter_if_three_or_more_last` separator.

    Parameters
    ----------
    texts : Tuple[str]
        List of all strings to be joined.
    delimiter_if_two : str
        String separating each element of `texts` if `len(texts) == 2`.
    delimiter_if_three_or_more_nonlast : str
        String separating each element _except_ the last two of `texts` if
        `len(texts) >= 3`.
    delimiter_if_three_or_more_last : str
        String separating the last two elements of `texts` if `len(texts) >= 3`.

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

# ....................{ JOINERS ~ as : conjunction         }....................
def join_as_conjunction(*texts) -> str:
    '''
    Conjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with `, `.
    * The last two passed strings are joined with `, and `.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' and ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', and '
    )


def join_as_conjunction_double_quoted(*texts) -> str:
    '''
    Conjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with `, `.
    * The last two passed strings are joined with `, and `.
    '''

    # Tuple of all passed strings double-quoted. Since the "*" operator applied
    # to this tuple below requires a sequence rather than generator, this is the
    # most space-efficient available sequence (i.e., frozen tuple).
    texts_quoted = tuple(double_quote(text) for text in texts)

    # Conjunctively join these strings.
    return join_as_conjunction(*texts_quoted)

# ....................{ JOINERS ~ as : conjunction         }....................
def join_as_disjunction(*texts) -> str:
    '''
    Disjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with `, `.
    * The last two passed strings are joined with `, or `.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' or ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', or '
    )


def join_as_disconjunction_double_quoted(*texts) -> str:
    '''
    Disjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with `, `.
    * The last two passed strings are joined with `, or `.
    '''

    # Tuple of all passed strings double-quoted. See
    # join_as_conjunction_double_quoted().
    texts_quoted = tuple(double_quote(text) for text in texts)

    # Disjunctively join these strings.
    return join_as_disjunction(*texts_quoted)

# ....................{ QUOTERS                            }....................
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

    . All prefixing and suffixing whitespace is removed from this string.
    . If this string is either double- or single-quoted, these quotes are
      removed.
    . All other double quotes in this string are escaped (e.g., replacing each
      `"` character with `\\"`).
    . The resulting string is double-quoted and returned.
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

# ....................{ REMOVERS                           }....................
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

# ....................{ REMOVERS ~ prefix                  }....................
@type_check
def remove_prefix(text: str, prefix: str, exception_message: str = None) -> str:
    '''
    Passed string with the passed prefix removed if present *or* raise an
    exception with the passed message (defaulting to a message synthesized from
    the passed arguments) otherwise.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    prefix : str
        Prefix to remove from this string.
    exception_message : optional[str]
        Exception message to be raised. Defaults to ``None``, in which case an
        exception message synthesized from the passed arguments is raised.

    Returns
    ----------
    str
        Resulting string as described above.
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
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    prefix : str
        Prefix to remove from this string.

    Returns
    ----------
    str
        Resulting string as described above.
    '''

    return text[len(prefix):] if is_prefix(text, prefix) else text

# ....................{ REMOVERS ~ prefix                  }....................
@type_check
def remove_suffix_if_found(text: str, suffix: str) -> str:
    '''
    Passed string with the passed suffix removed if present _or_ the passed
    string as is otherwise.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    suffix : str
        Suffix to remove from this string.

    Returns
    ----------
    str
        Resulting string as described above.
    '''

    # There exists a special case *NOT* present in remove_prefix(). If this
    # suffix is empty, "string[:-0]" is also incorrectly empty. Avoid returning
    # the empty string in this case by explicitly testing for emptiness.
    return text[:-len(suffix)] if suffix and is_suffix(text, suffix) else text


@type_check
def remove_suffix_with_prefix(text: str, suffix_prefix: str) -> str:
    '''
    Return the passed string without the first instance of the passed suffix
    prefix and all characters following this suffix prefix in this string if
    present _or_ raise an exception otherwise.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    suffix_prefix : str
        Non-empty substring of this string to begin removing characters at.

    Returns
    ----------
    str
        Resulting string as described above.

    Raises
    ----------
    BetseStrException
        If either:
        * This string does _not_ contain this suffix prefix.
        * This suffix prefix is the empty string.
    '''

    # If this suffix prefix is the empty string, raise an exception.
    die_if_empty(text=text, exception_message='Suffix prefix empty.')

    # Attempt to return the desired string.
    try:
        return text[:text.index(suffix_prefix)]
    # If doing so fails with a builtin exception...
    except ValueError as exception:
        # ...failing to provide the contents of these arguments, wrap this
        # exception with a fine-grained exception providing these contents.
        if str(exception) == 'substring not found':
            raise BetseStrException(
                'Text "{}" contains no suffix prefix "{}".'.format(
                    text, suffix_prefix))

        # Else, re-raise this exception as is.
        raise

# ....................{ CASERS                             }....................
#FIXME: For orthogonality, rename to uppercase_char_first().
@type_check
def uppercase_first_char(text: str) -> str:
    '''
    Uppercase the first character of the passed string.

    Whereas the related :meth:`str.capitalize` method both uppercases the first
    character of this string *and* lowercases all remaining characters, this
    function *only* uppercases the first character. All remaining characters
    remain unmodified.
    '''

    return (
        text[0].upper() + (text[1:] if len(text) > 2 else '')
        if len(text) else ''
    )

# ....................{ QUOTERS                            }....................
@type_check
def shell_quote(text: str) -> str:
    '''
    Shell-quote the passed string in a platform-specific manner.

    If the current platform is:

    * _Not_ Windows (e.g., Linux, OS X), the returned string is guaranteed to be
      suitable for passing as an arbitrary positional argument to external
      commands.
    * Windows, the returned string is suitable for passing _only_ to external
      commands parsing arguments in the same manner as the Microsoft C runtime.
      While _all_ applications on POSIX-compliant systems are required to parse
      arguments in the same manner (i.e., according to Bourne shell lexing), no
      such standard applies to Windows applications. Shell quoting is therefore
      fragile under Windows -- like pretty much everything.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current OS is Windows, do *NOT* perform POSIX-compatible quoting.
    # Windows is POSIX-incompatible and hence does *NOT* parse command-line
    # arguments according to POSIX standards. In particular, Windows does *NOT*
    # treat single-quoted arguments as single arguments but rather as multiple
    # shell words delimited by the raw literal `'`. This is circumventable by
    # calling an officially undocumented Windows-specific function. (Awesome.)
    if oses.is_windows():
        import subprocess
        return subprocess.list2cmdline([text])
    # Else, perform POSIX-compatible quoting.
    else:
        import shlex
        return shlex.quote(text)

# ....................{ WRAPPERS                           }....................
def wrap_lines(lines: IterableTypes, **kwargs) -> str:
    '''
    Wrap the passed iterable of lines to the passed line width, prefixing each
    resulting wrapped line by the passed line prefix.

    See Also
    ----------
    wrap()
        For further details.
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

    * ``text_wrapper`, the object on which to call the ``wrap()`` function or
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
        'Text wrapper "{}" attribute "wrap" not callable.'.format(text_wrapper))

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

# ....................{ WRAPPERS ~ un                      }....................
@type_check
def unwrap(text: str) -> str:
    '''
    Convert the passed possibly **multiline string** (i.e., string containing
    one or more newlines) to a **single-line string** (i.e., string containing
    no newlines).
    '''

    # Nice one, stdlib.
    return text.replace('\n', ' ')

# ....................{ (IN|DE)DENTERS                     }....................
def dedent(*texts) -> str:
    '''
    Remove all indentation shared in common by all lines of all passed strings.
    '''

    return textwrap.dedent(*texts)
