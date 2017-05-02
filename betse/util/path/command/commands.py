#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **command** (i.e., external executable file) facilities.
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.exceptions import BetseCommandException
from betse.metadata import SCRIPT_BASENAME
from betse.util.type.types import type_check

# ....................{ GLOBALS                            }....................
_CURRENT_BASENAME = None
'''
Basename of the command originating the active Python interpreter.

This global caches the return value of the frequently called
:func:`get_current_basename` function, for efficiency.
'''

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_command(pathname: str) -> None:
    '''
    Raise an exception unless a command with the passed path exists.

    Parameters
    ----------
    pathname : str
        Basename or absolute or relative path of the executable file to inspect.

    See Also
    ----------
    :func:`is_command`
        Further details.
    '''

    if not is_command(pathname):
        raise BetseCommandException(
            '"{}" not an executable command.'.format(pathname))

# ....................{ TESTERS                            }....................
@type_check
def is_command(pathname: str) -> None:
    '''
    `True` only if a command with the passed path exists.

    This is the case if this path is either:

    * The basename of an executable file in the current `${PATH}`.
    * The relative or absolute path of an executable file.

    Parameters
    ----------
    pathname : str
        Basename or absolute or relative path of the executable file to inspect.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import files, paths
    from betse.util.path.command import pathables

    # This path is that of an existing command if and only if either...
    return (
        # This path is that of an executable file *OR*
        files.is_file_executable(pathname) or (
        # This path is that of a basename in the current ${PATH}.
        paths.is_basename(pathname) and pathables.is_pathable(pathname))
    )

# ....................{ GETTERS                            }....................
def get_current_basename() -> str:
    '''
    Get the basename of the command originating the active Python interpreter.

    If this interpreter is interpreting a block of arbitrary runtime code passed
    to this interpreter on the command line via Python's `-c` option (e.g., due
    to being called by a distributed `pytest-xdist` test), this function
    unconditionally returns the basename of the BETSE CLI (e.g., `betse`) rather
    than `-c`. `-c` is a rather unexpected and non-human-readable basename!
    '''

    # If this function has already been called, return the string cached by the
    # first call to this function.
    global _CURRENT_BASENAME
    if _CURRENT_BASENAME is not None:
        return _CURRENT_BASENAME

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Raw absolute or relative path of the current command.
    _CURRENT_BASENAME = sys.argv[0]

    # If this is the non-human-readable "-c" Python interpreter option,
    # substitute this with the human-readable basename of the BETSE CLI.
    if _CURRENT_BASENAME == '-c':
        _CURRENT_BASENAME = SCRIPT_BASENAME
    # Else, reduce this absolute or relative path to a basename.
    else:
        _CURRENT_BASENAME = paths.get_basename(_CURRENT_BASENAME)

    # Return this basename.
    return _CURRENT_BASENAME
