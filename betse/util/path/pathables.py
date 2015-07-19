#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level *pathable* (i.e., external commands in the current `${PATH}`)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionPath
from betse.util.path import paths
import shutil

# ....................{ TESTERS                            }....................
def is_pathable(command_basename: str) -> bool:
    '''
    True if the external command with the passed basename exists.

    Specifically, return True if such basename is that of an executable file in
    the current `${PATH}`. If such basename contains a directory separator and
    is hence *not* a basename, an exception is raised.
    '''
    assert isinstance(command_basename, str),\
        '"{}" not a string.'.format(command_basename)
    assert len(command_basename), 'Command name empty.'

    # If such string is *NOT* a pure basename, fail.
    if paths.is_dirname_empty(command_basename):
        raise BetseExceptionPath(
            'Command "{}" contains directory separators.'.format(
                command_basename))

    # Return whether such command is found.
    return shutil.which(command_basename) is not None

# --------------------( WASTELANDS                         )--------------------
