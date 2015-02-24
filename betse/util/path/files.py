#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory file facilities.

This module is named `files` rather than `file` to avoid conflict with the stock
`files` class.
'''

# ....................{ IMPORTS                            }....................
from os import path
import os

from betse.exceptions import BetseExceptionFile


# ....................{ EXCEPTIONS                         }....................
def die_unless_found(filename: str) -> None:
    '''
    Raise an exception unless the passed non-directory file exists.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    if not is_file(filename):
        raise BetseExceptionFile(
            'File "{}" not found or not a readable file.'.format(filename))

# ....................{ TESTERS                            }....................
def is_file(filename: str) -> bool:
    '''
    True if the passed non-directory file exists.

    Versus path.isfile()
    ----------
    This function fundamentally differs from the stock `path.isfile()` function.
    Whereas the latter returns True only for non-special files and hence False
    for all non-directory special files (e.g., device nodes, symbolic links),
    this function returns True for for *all* non-directory files regardless of
    whether such files are special or not.

    Why? Because this function complies with POSIX semantics, whereas
    `path.isfile()` does *not*. Under POSIX, it is largely irrelevant whether a non-directory
    file is special or not; it simply matters whether such file is a directory
    or not. For example, the external command `rm` removes only non-directory
    files and the external command `rmdir` removes only empty directories.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    return path.exists(filename) and not path.isdir(filename)

# ....................{ REMOVERS                           }....................
def remove(filename: str) -> None:
    '''
    Remove the passed file.

    If such file does *not* exist, an exception is raised.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    os.remove(filename)

# --------------------( WASTELANDS                         )--------------------
