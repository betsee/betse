#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
**Command** (i.e., external executable file) facilities.
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.metadata import SCRIPT_NAME_CLI

# ....................{ GLOBALS                            }....................
_CURRENT_BASENAME = None
'''
Basename of the command originating the active Python interpreter.

This global caches the return value of the frequently called
`get_current_basename()` function, for efficiency.
'''

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
        _CURRENT_BASENAME = SCRIPT_NAME_CLI
    # Else, reduce this absolute or relative path to a basename.
    else:
        _CURRENT_BASENAME = paths.get_basename(_CURRENT_BASENAME)

    # Return this basename.
    return _CURRENT_BASENAME
