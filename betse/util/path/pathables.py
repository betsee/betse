#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level *pathable* (i.e., external commands in the current `${PATH}`)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionPath
from betse.util.type import types
import shutil

# ....................{ TESTERS                            }....................
def is_pathable(command_basename: str) -> bool:
    '''
    `True` if the external command with the passed basename exists _or_ `False`
    otherwise.

    Specifically, `True` if this basename is that of an executable file in the
    current `${PATH}`. If this basename contains a directory separator and is
    hence _not_ a basename, an exception is raised.
    '''
    assert types.is_str_nonempty(command_basename), (
        types.assert_not_str_nonempty(command_basename, 'Command name'))

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # If this string is *NOT* a pure basename, fail.
    if paths.is_basename(command_basename):
        raise BetseExceptionPath(
            'Command "{}" contains directory separators.'.format(
                command_basename))

    # Get whether this command exists or not.
    return shutil.which(command_basename) is not None

# --------------------( WASTELANDS                         )--------------------
