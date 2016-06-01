#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

# ....................{ IMPORTS                            }....................
import textwrap
from betse.exceptions import BetseExceptionString
from betse.util.type import types
from textwrap import TextWrapper

# ....................{ SINGLETONS                         }....................
text_wrapper = TextWrapper()
'''
Singleton `TextWrapper` instance with which to wrap text.

Such singleton improves efficiency -- occasionally dramatically. *All* public
functions provided by module `textwrap` implicitly instantiate temporary
`TextWrapper` instances on each call to such functions.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_prefix(
    text: str, prefix: str, exception_message: str = None) -> None:
    '''
    Raise an exception with the passed message (defaulting to a message
    synthesized from the passed arguments) if the second passed string does
    _not_ prefix the first passed string.

    Parameters
    ----------
    text : str
        String to be tested.
    prefix : str
        Prefix to test for.
    exception_message : optional[str]
        Exception message to be raised. Defaults to `None`, in which case an
        exception message synthesized from the passed arguments is raised.
    '''

    # If this string is *NOT* prefixed by this prefix, raise an exception.
    if not is_prefix(text, prefix):
        # If no exception message was passed, synthesize one from this name.
        if not exception_message:
            exception_message = 'Text "{}" not prefixed by "{}".'.format(
                text, prefix)
        assert types.is_str(exception_message), (
            types.assert_not_str(exception_message))

        # Raise this exception.
        raise BetseExceptionString(exception_message)

# ....................{ TESTERS                            }....................
def is_prefix(text: str, prefix: str) -> bool:
    '''
    `True` only if the second passed string prefixes the first passed string.

    Parameters
    ----------
    text : str
        String to be tested.
    prefix : str
        Prefix to test for.

    Returns
    ----------
    bool
        `True` only if `prefix` prefixes `text`.
    '''
    assert types.is_str(text), types.assert_not_str(text)
    assert types.is_str_nonempty(prefix), (
        types.assert_not_str_nonempty(prefix, 'Prefix'))

    return text.startswith(prefix)


def is_suffix(text: str, suffix: str) -> bool:
    '''
    `True` only if the second passed string suffixes the first passed string.

    Parameters
    ----------
    text : str
        String to be tested.
    suffix : str
        Prefix to test for.

    Returns
    ----------
    bool
        `True` only if `suffix` suffixes `text`.
    '''
    assert types.is_str(text), types.assert_not_str(text)
    assert types.is_str_nonempty(suffix), (
        types.assert_not_str_nonempty(suffix, 'Suffix'))

    return text.endswith(suffix)

# ....................{ JOINERS                            }....................
def join(*texts) -> str:
    '''
    Concatenate the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard `"".join((...))`
    method, whose syntax is arguably overly verbose.
    '''
    return join_on(*texts, delimiter='')


def join_on_newline(*texts) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    `"\n".join((...))` method, whose syntax is arguably overly verbose.
    '''
    return join_on(*texts, delimiter='\n')


def join_on(*texts, delimiter: str) -> str:
    '''
    Join the passed strings with the passed separating delimiter.

    This is a convenience function wrapping the standard
    `"...".join((...))` method, whose syntax is arguably overly obfuscated.
    '''
    # To avoid obscure chicken-and-egg exceptions when logging exceptions, this
    # delimiter must be validated manually rather than calling types.is_str().
    assert isinstance(delimiter, str), types.assert_not_str(delimiter)

    # If only one object was passed and such object is a non-string iterable
    # (e.g, list, tuple), set the list of passed strings to such object.
    if len(texts) == 1 and types.is_iterable_nonstr(texts[0]):
        texts = texts[0]

    # Join such texts.
    return delimiter.join(texts)

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

# ....................{ REMOVERS                           }....................
def remove_presuffix_whitespace(text: str) -> str:
    '''
    New string with all prefixing and suffixing whitespace removed from the
    passed string.

    If this string is neither prefixed nor suffixed by whitespace, this string
    is returned as is.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.

    Returns
    ----------
    str
        Resulting string as described above.
    '''
    assert types.is_str(text), types.assert_not_str(text)

    return text.strip()

# ....................{ REMOVERS ~ prefix                  }....................
def remove_prefix(text: str, prefix: str, exception_message: str = None) -> str:
    '''
    New string with the passed prefix removed from the passed string if the
    former prefixes the latter _or_ raise an exception with the passed message
    (defaulting to a message synthesized from the passed arguments) otherwise.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    prefix : str
        Prefix to remove from this string.
    exception_message : optional[str]
        Exception message to be raised. Defaults to `None`, in which case an
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


def remove_prefix_if_found(text: str, prefix: str) -> str:
    '''
    New string with the passed prefix removed from the passed string if the
    former prefixes the latter _or_ the passed string as is otherwise.

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
def remove_suffix_if_found(text: str, suffix: str) -> str:
    '''
    New string with the passed suffix removed from the passed string if the
    former suffixes the latter _or_ the passed string as is otherwise.

    Parameters
    ----------
    text : str
        String to be examined. Since strings are immutable in Python, this
        string remains unmodified.
    suffix : str
        Prefix to remove from this string.

    Returns
    ----------
    str
        Resulting string as described above.
    '''

    # There exists a special case *NOT* present in remove_prefix(). If this
    # suffix is empty, "string[:-0]" is also incorrectly empty. Avoid returning
    # the empty string in this case by explicitly testing for emptiness.
    return text[:-len(suffix)] if suffix and is_suffix(text, suffix) else text

# ....................{ CASERS                             }....................
def uppercase_first_char(text: str) -> str:
    '''
    Uppercase the first character of the passed string.

    Whereas the related `str.capitalize()` method both uppercases the first
    character of this string _and_ lowercases all remaining characters, this
    function _only_ uppercases the first character. All remaining characters
    remain unmodified.
    '''
    assert types.is_str(text), types.assert_not_str(text)

    return (
        text[0].upper() + (text[1:] if len(text) > 2 else '')
        if len(text) else ''
    )

# ....................{ QUOTERS                            }....................
def shell_quote(text: str) -> str:
    '''
    Shell-quote the passed string.

    If the current operating system is:

    * *Not* Windows (e.g., Linux, OS X), the returned string is guaranteed to be
      suitable for passing as an arbitrary positional argument to external
      commands.
    * Windows, the returned string is suitable for passing *only* to external
      commands parsing arguments according in the same manner as the Microsoft C
      runtime. Whereas *all* applications running under POSIX-compliant systems
      are required to parse arguments in the same manner (e.g., according to
      Bourne shell lexing), no such standard applies to applications running
      under Windows. For this reason, shell quoting is inherently unreliable
      under Windows.
    '''
    assert types.is_str(text), types.assert_not_str(text)

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current OS is Windows, do *NOT* perform POSIX-compatible quoting.
    # Windows is POSIX-incompatible and hence does *NOT* parse command-line
    # arguments according to POSIX standards. In particular, Windows does *NOT*
    # treat single-quoted arguments as single arguments but rather as multiple
    # shell words delimited by the raw literal `'`.  This is circumventable by
    # calling an officially undocumented Windows-specific function. (Awesome.)
    if oses.is_windows():
        import subprocess
        return subprocess.list2cmdline([text])
    # Else, perform POSIX-compatible quoting.
    else:
        import shlex
        return shlex.quote(text)

# ....................{ WRAPPERS                           }....................
def wrap_lines(lines: list, **kwargs) -> str:
    '''
    Wrap the passed iterable of lines to the passed line width, prefixing each
    resulting wrapped line by the passed line prefix.

    See Also
    ----------
    wrap()
        For further details.
    '''
    return wrap(join(lines), **kwargs)


def wrap(
    text: str,
    text_wrapper = textwrap,
    line_prefix: str = '',
    **kwargs) -> str:
    '''
    Wrap the passed text to the passed line width, prefixing each resulting
    wrapped line by the passed line prefix.

    This function accepts the following keyword arguments in addition to those
    accepted by `textwrap.wrap()`:

    * `text_wrapper`, the object on which to call the `wrap()` function or
      method. For safety, this defaults to module `textwrap`.
    * `line_prefix`, the substring prefixing each wrapped output line.

    See Also
    ----------
    https://docs.python.org/3/library/textwrap.html
        For further details on keyword arguments.
    '''
    assert types.is_str(text), types.assert_not_str(text)
    assert types.is_str(line_prefix), types.assert_not_str(line_prefix)
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

# ....................{ (IN|DE)DENTERS                     }....................
def dedent(*texts) -> str:
    '''
    Remove all indentation shared in common by all lines of all passed strings.
    '''
    return textwrap.dedent(*texts)
