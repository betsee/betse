#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level temporary path facilities, including both temporary directories _and_
non-directory files.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, StrOrNoneTypes
from tempfile import NamedTemporaryFile

# ....................{ WRITERS                            }....................
# The type of the return value is the private class
# "tempfile._TemporaryFileWrapper", which due to being private is intentionally
# *NOT* type-checked here.
@type_check
def write_bytes(encoding: StrOrNoneTypes = None):
    '''
    Open and return a temporary named binary file for byte-oriented writing.

    Parameters
    ----------
    encoding : optional[str]
        Name of the encoding with which to encode bytes into plaintext strings
        if desired or `None` otherwise (i.e., if bytes are to be written as is
        without being encoded into such strings). Defaults to `None`.

    Returns
    ----------
    tempfile._TemporaryFileWrapper
        `file`-like object encapsulating the opened file.

    See Also
    ----------
    :func:`writing_chars`
        Further details.
    '''

    return NamedTemporaryFile(delete=False, encoding=encoding)


# The type of the return value is the private class
# "tempfile._TemporaryFileWrapper", which due to being private is intentionally
# *NOT* type-checked here.
@type_check
def write_chars(encoding: str = 'utf-8'):
    '''
    Open and return a temporary named plaintext file for line-oriented writing.

    This function returns a `file`-like object, suitable for use wherever the
    builtin `open()` would otherwise be called (e.g., in `with` statements).
    The filename for the underlying file is accessible in the typical way for
    named `file`-like objects (namely, via the `name` attribute of such object).
    The underlying file will _not_ be implicitly deleted when such object is
    closed, which typically defeats the purpose of creating a temporary named
    file in the first place.

    Text Lines
    ----------
    Temporary files are _not_ reliably openable for line-oriented writing in a
    cross-platform manner. Hence, callers are discouraged from calling this
    function except where absolutely required (e.g., some or all underlying
    encodings are unknown). Instead, text lines may be manually written in a
    byte-oriented manner by encoding such lines under a known encoding _and_
    appending such lines by the newline delimiter specific to the current
    platform: e.g.,

        >>> from betse.util.path imports temps
        >>> import os
        >>> tempfile = temps.writing_bytes(encoding='utf-8')
        >>> tempfile.write(b'The strength of a civilization is not measured ')
        >>> tempfile.write(b'by its ability to fight wars, but rather by its ')
        >>> tempfile.write(b'ability to prevent them.')
        >>> tempfile.write(bytes(os.linesep, encoding='utf-8'))

    Parameters
    ----------
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    tempfile._TemporaryFileWrapper
        `file`-like object encapsulating the opened file.
    '''

    return NamedTemporaryFile(mode='w+', delete=False, encoding=encoding)
