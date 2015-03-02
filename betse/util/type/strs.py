#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

# ....................{ IMPORTS                            }....................
from textwrap import TextWrapper
import textwrap

# ....................{ SINGLETONS                         }....................
text_wrapper = TextWrapper()
'''
Singleton `TextWrapper` instance with which to wrap text.

Such singleton improves efficiency -- occasionally dramatically. *All* public
functions provided by module `textwrap` implicitly instantiate temporary
`TextWrapper` instances on each call to such functions.
'''

# ....................{ TESTERS                            }....................
def is_str(obj) -> bool:
    '''
    True if the passed object is a *string* (i.e., is an instance of the `str`
    class or a subclass thereof).
    '''
    return isinstance(obj, str)

# ....................{ JOINERS                            }....................
def join(*texts) -> str:
    '''
    Concatenate the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard `"".join((...))`
    method, whose syntax is arguably overly verbose.
    '''
    return join_on(*texts, delimiter = '')

def join_on_newline(*texts) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    `"\n".join((...))` method, whose syntax is arguably overly verbose.
    '''
    return join_on(*texts, delimiter = '\n')

def join_on(*texts, delimiter: str) -> str:
    '''
    Join the passed strings with the passed separating delimiter.

    This is a convenience function wrapping the standard
    `"...".join((...))` method, whose syntax is arguably overly verbose.
    '''
    assert isinstance(delimiter, str), '"{}" not a string.'.format(delimiter)

    # If only one object was passed and...
    if len(texts) == 1:
        # Avoid circular import dependencies.
        from betse.util.type import containers

        # ...such object is a non-string iterable (e.g, list, tuple), set the
        # list of passed strings to such object.
        if containers.is_iterable_nonstring(texts[0]):
            texts = texts[0]

    # Join such texts.
    return delimiter.join(texts)

# ....................{ REMOVERS                           }....................
def remove_prefix_if_found(text: str, prefix: str) -> str:
    '''
    Remove the passed prefix from the passed string if such string is prefixed
    by such prefix *or* such string as is otherwise.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    assert isinstance(prefix, str), '"{}" not a string.'.format(prefix)
    return text[len(prefix):] if text.startswith(prefix) else text

def remove_suffix_if_found(text: str, suffix: str) -> str:
    '''
    Remove the passed suffix from the passed string if such string is suffixed
    by such suffix *or* such string as is otherwise.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    assert isinstance(suffix, str), '"{}" not a string.'.format(suffix)

    # There exists a special case *NOT* present in remove_prefix(). If such
    # suffix is empty, "string[:-0]" is also incorrectly empty. Avoid returning
    # the empty string in such case by explicitly testing for emptiness.
    return text[:-len(suffix)] if suffix and text.endswith(suffix) else text

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

def wrap(text: str, line_prefix: str = '', **kwargs) -> str:
    '''
    Wrap the passed text to the passed line width, prefixing each resulting
    wrapped line by the passed line prefix.

    This function accepts the same keyword arguments as `textwrap.wrap()`, in
    addition to the `line_prefix` argument specific to this function.

    See Also
    ----------
    https://docs.python.org/3/library/textwrap.html
        For further details on keyword arguments.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)
    assert isinstance(line_prefix, str),\
        '"{}" not a string.'.format(line_prefix)

    # If passed a nonempty line prefix, add appropriate keyword arguments.
    if line_prefix:
        kwargs['initial_indent'] = line_prefix
        kwargs['subsequent_indent'] = line_prefix

    # For reliability, call the wrap() function of module "textwrap" rather than
    # the singleton object "text_wrapper". The former internally instantiates a
    # new instance of class "TextWrapper" and hence resets all wrapping
    # attributes to sensible defaults, whereas the latter reuses existing such
    # attributes -- which may no longer retain sensible defaults.
    return join_on_newline(textwrap.wrap(text, **kwargs))

# ....................{ (IN|DE)DENTERS                     }....................
def dedent(*texts) -> str:
    '''
    Remove all indentation shared in common by all lines of all passed strings.
    '''
    return textwrap.dedent(*texts)

# --------------------( WASTELANDS                         )--------------------
# (defaulting to the empty string)
