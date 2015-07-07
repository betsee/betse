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
from betse.exceptions import BetseExceptionPath
from os import path

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_dirname_empty(
    pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed pathname is a pure basename.

    See Also
    ----------
    `is_dirname_empty()`
        For further details.
    '''
    if not is_dirname_empty(pathname):
        # If no message was passed, default such message.
        if not exception_message:
            exception_message =\
                'Path "{}" contains directory separators.'.format(pathname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise BetseExceptionPath(exception_message)

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_special(pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception if the passed path is an existing special path.

    See Also
    ----------
    `is_special()`
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
        raise BetseExceptionPath(exception_message)

# ....................{ TESTERS                            }....................
def is_path(pathname: str) -> bool:
    '''
    True if the passed path exists.

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
    True if the passed path is an existing special file.

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

# ....................{ TESTERS ~ pathname                 }....................
def is_dirname_empty(pathname: str) -> bool:
    '''
    True if the passed pathname is a *pure basename* (i.e., contains no
    directory separators and hence no directory components).
    '''
    return path.sep in pathname

def is_filetype(pathname: str, filetype: str) -> bool:
    '''
    True if the passed pathname has the passed filetype.

    Such filetype may contain arbitrarily many `.` characters, including an
    optional prefixing `.`. Regardless, this function behaves as expected.
    '''
    assert isinstance(filetype, str), '"{}" not a string.'.format(filetype)

    # Avoid circular import dependencies.
    from betse.util.type import strs

    # Test such filetype, prefixed by "." unless already prefixed.
    return strs.is_suffix(
        pathname,
        strs.add_prefix_unless_found(filetype, '.'))

# ....................{ GETTERS                            }....................
def get_basename(pathname: str) -> str:
    '''
    Get the **basename** (i.e., last component) of the passed path.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.basename(pathname)

def get_filetype(pathname: str) -> str:
    '''
    Get the last **filetype** (i.e., last `.`-prefixed substring of the
    basename *not* including such `.`) of the passed path if such path has a
    filetype or None otherwise.

    If such path has multiple filetypes (e.g., `odium.reigns.tar.gz`), only the
    last such filetype is returned.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Such filetype. (Yes, splitext() is exceedingly poorly named.)
    filetype = path.splitext(pathname)[1]

    # Get such filetype, stripping the prefixing "." from the string returned by
    # the prior call if such path has a filetype or otherwise returning None.
    return filetype[1:] if filetype else None

# ....................{ GETTERS ~ dirname                  }....................
def get_dirname(pathname: str) -> str:
    '''
    Get the *dirname* (i.e., parent directory) of the passed path if such path
    has a dirname or raise an exception otherwise.
    '''
    # Ensure such path has a dirname.
    die_unless_dirname_empty(pathname)

    # Get such dirname. Technically, the above call *SHOULD* have ensured such
    # dirname to exist. This is a sufficiently critical function, however, to
    # warrant asserting this constraint for safety.
    dirname = get_dirname_or_empty(pathname)
    assert len(dirname), 'Pathname "{}" dirname empty.'.format(pathname)
    return dirname

def get_dirname_or_empty(pathname: str) -> str:
    '''
    Get the *dirname* (i.e., parent directory) of the passed path if such path
    has a dirname or the empty string otherwise.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.dirname(pathname)

# ....................{ REMOVERS                           }....................
def remove_filetype_if_found(pathname: str) -> str:
    '''
    Remove the last filetype (including prefixing `.`) from the passed path if
    such path has a filetype *or* return such path as is otherwise.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.splitext(pathname)[0]

# ....................{ JOINERS                            }....................
def join(*pathnames) -> str:
    '''
    Join the passed pathnames on the directory separator specific to the current
    operating system.

    This is a convenience function wrapping the standard `os.path.join()`
    function _without_ adding functionality to such function -- principally to
    unify and hence simplify `import` statements in other modules.
    '''
    return path.join(*pathnames)

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

# --------------------( WASTELANDS                         )--------------------
# def get_dirname_or_none(pathname: str) -> str:
#     '''
#     Get the *dirname* (i.e., parent directory) of the passed path if such path
#     has a dirname or None otherwise.
#     '''
#     return get_dirname_or_empty(pathname) or None

    # Strip the prefixing "." from such filetype if any.
    # filetype = strs.remove_prefix_if_any(filetype, '.')
    #
    # # Test such filetype.
    # return get_filetype(pathname) == filetype

    # Such filetype may be either prefixed by `.` *or* not prefixed by `.`, but
    # should otherwise contain *no* `.` characters (e.g., `.gz` is accepted but
    # `.tar.gz` is *not*). In either case, this function operates as expected.

    # While there are more efficient implementations, the simplest should be
    # fine... for now.
    # return pathname == path.get_basename(pathname)

    # assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    # assert len(pathname), 'Pathname empty.'
    #
    # # Such dirname.
    # dirname = path.dirname(pathname)
    #
    # # Get such dirname. Since path.dirname() returns the empty string rather
    # # than None for paths without a dirname, convert the former to the latter.
    # return dirname if dirname else None
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
