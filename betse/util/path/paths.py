#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level path facilities specific to neither directories nor non-directory
files.

This module is named `paths` rather than `path` to avoid conflict with the stock
`path` module of the `os` package.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionFile
from os import path

# ....................{ EXCEPTIONS                         }....................
def die_if_special(pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception if the passed path is an existing special path.

    See Also
    ----------
    is_special
        For further details.
    '''
    if is_special(pathname):
        # If no message was passed, default such message.
        if not exception_message:
            if path.isdir(pathname):
                exception_message =\
                    'Path "{}" already an existing directory.'.format(pathname)
            elif path.islink(pathname):
                exception_message =\
                    'Path "{}" already an existing symbolic link.'.format(
                        pathname)
            else:
                exception_message = 'Path "{}" already a special file.'.format(
                    pathname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise BetseExceptionFile(exception_message)

# ....................{ TESTERS                            }....................
def is_path(pathname: str) -> bool:
    '''
    True if the passed exists.

    If such path is an existing **broken symbolic link** (i.e., a symbolic link
    whose target no longer exists), this function still returns True.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Since path.exists() returns False for broken symbolic links, defer to
    # path.lexists() instead.
    return path.lexists(pathname)

def is_special(pathname: str) -> bool:
    '''
    True if the passed is an existing special file.

    Special files include directories, device nodes, sockets, and symbolic
    links.
    '''
    # True if such path exists and...
    return is_path(pathname) and (
        # ...is either a symbolic link *OR* neither a regular file nor symbolic
        # link to such a file. In the latter case, predicate logic guarantees
        # such file to *NOT* be a symbolic link, thus reducing this test to:
        # "...is either a symbolic link *OR* not a regular file."
        path.islink(pathname) or\
        not path.isfile(pathname)
    )

# ....................{ GETTERS                            }....................
def get_dirname(pathname: str) -> str:
    '''
    Get the *dirname* (i.e., parent directory) of the passed path if such path
    has a dirname or None otherwise.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Such dirname.
    dirname = path.dirname(pathname)

    # Get such dirname. Since path.dirname() returns the empty string rather
    # than None for paths without a dirname, convert the former to the latter.
    return dirname if dirname else None

def get_basename(pathname: str) -> str:
    '''
    Get the **basename** (i.e., last component) of the passed path.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.basename(pathname)

def get_filetype(pathname: str) -> str:
    '''
    Get the **filetype** (i.e., last `.`-prefixed substring of the basename) of
    the passed path has a filetype or None otherwise.

    If such has multiple filetypes (e.g., `odium.reigns.tar.gz`), only the last
    such filetype is returned.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Such filetype. (Yes, splitext() is exceedingly poorly named.)
    filetype = path.splitext(pathname)[1]

    # Get such filetype, stripping the prefixing "." from the string returned by
    # the prior call if such path has a filetype or otherwise returning None.
    return filetype[1:] if filetype else None

# ....................{ CANONICALIZERS                     }....................
def canonicalize(pathname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed path.

    Specifically (in order):

    * Perform **tilde expansion,** replacing a `~` character prefixing such path
      by the absolute path of the current user's home directory.
    * Perform **path normalization,** thus:
      * Collapsing redundant separators (e.g., converting `//` to `/`).
      * Converting relative to absolute path components (e.g., converting `../`
        to the name of the parent directory of such component).
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.abspath(path.expanduser(pathname))

# ....................{ JOINERS                            }....................
def join(*pathnames) -> str:
    '''
    Join the passed pathnames on the directory separator specific to the current
    operating system.

    This is a convenience function wrapping the standard `os.path.join()`
    function, provided to reduce the number of import statements required by
    other modules.
    '''
    return path.join(*pathnames)

# --------------------( WASTELANDS                         )--------------------
#FUXME: Obviously insufficient. This should also test whether such path is a

# ....................{ REMOVERS                           }....................
# def remove(filename: str) -> None:
#     '''
#     Remove the passed file.
#
#     If such file does *not* exist, an exception is raised.
#     '''
#     assert isinstance(filename, str), '"{}" not a string.'.format(filename)
#     os.remove(filename)

# def __init__():
#     '''
#
#     To support caller-specific exception handling, this function does *not*
#     validate such constants (e.g., by creating non-existent directories). See
#     `paths.init()` for such functionality.
#     '''
#     # Declare such constants to be globals, permitting their modification below.
#     global DATA_DIRNAME, DOT_DIRNAME
