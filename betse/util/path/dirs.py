#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level directory facilities.

This module is named `dirs` rather than `dir` to avoid conflict with the `dir()`
builtin.
'''

# ....................{ IMPORTS                            }....................
from betse.util.exceptions import BetseExceptionDir
from os import path
import os

# ....................{ EXCEPTIONS                         }....................
def die_unless_found(dirname: str) -> None:
    '''
    Raise an exception unless the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    if not is_dir(dirname):
        raise BetseExceptionDir(
            'Directory "{}" not found or not a readable directory.'.format(
                dirname))

def die_unless_parent_found(pathname: str) -> None:
    '''
    Raise an exception unless the parent directory of the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    die_unless_found(get_dirname(pathname))

# ....................{ TESTERS                            }....................
def is_dir(dirname: str) -> bool:
    '''
    True if the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    return path.isdir(dirname)

# ....................{ GETTERS                            }....................
def get_dirname(pathname: str) -> str:
    '''
    Get the *dirname* (i.e., parent directory) of the passed path.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.dirname(pathname)

# ....................{ MAKERS                             }....................
#FIXME: Replace all existing calls to os.makedirs() by calls to such functions.

def make_unless_found(dirname: str) -> None:
    '''
    Create the passed directory if such directory does *not* already exist.

    All nonexistent parents of such directory will also be recursively created,
    mimicking the action of the conventional shell command `mkdir -p`.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    os.makedirs(dirname, exist_ok = True)

def make_parent_unless_found(*pathnames) -> None:
    '''
    Create the parent directory of each passed path for parent directories that
    do *not* already exist.
    '''
    for pathname in pathnames:
        make_unless_found(get_dirname(pathname))

# --------------------( WASTELANDS                         )--------------------
# from betse.util.path import paths
