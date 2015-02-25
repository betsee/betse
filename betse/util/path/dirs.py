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
from betse.exceptions import BetseExceptionDir
from betse.util.io import loggers
from betse.util.path import paths
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
    die_unless_found(paths.get_dirname(pathname))

# ....................{ TESTERS                            }....................
def is_dir(dirname: str) -> bool:
    '''
    True if the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    return path.isdir(dirname)

# ....................{ MAKERS                             }....................
def make_unless_found(dirname: str) -> None:
    '''
    Create the passed directory if such directory does *not* already exist.

    All nonexistent parents of such directory will also be recursively created,
    mimicking the action of the conventional shell command `mkdir -p`.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    assert len(dirname), 'Dirname empty.'

    # If such directory does *NOT* already exist, create such directory. To
    # support logging, such condition is explicitly tested for. To avoid race
    # conditions (e.g., in the event such directory is created between testing
    # and creating such directory), we preserve the makedirs() keyword argument
    # "exist_ok = True".
    if not is_dir(dirname):
        # Log such creation.
        loggers.log_info('Creating directory "%s".', dirname)

        # Create such directory if still needed.
        os.makedirs(dirname, exist_ok = True)

def make_parent_unless_found(*pathnames) -> None:
    '''
    Create the parent directory of each passed path if such directory does *not*
    already exist.

    See Also
    ----------
    make_unless_found()
        For further details.
    '''
    # Canonicalize each pathname *BEFORE* attempting to get its dirname.
    # Relative pathnames do *NOT* have sane dirnames (e.g., the dirname for a
    # relative pathname "metatron" is the empty string) and hence *MUST* be
    # converted to absolute pathnames first.
    for pathname in pathnames:
        make_unless_found(
            paths.get_dirname(
                paths.canonicalize(pathname)))

# --------------------( WASTELANDS                         )--------------------
#FUXME: Replace all existing calls to os.makedirs() by calls to such functions.
# from betse.util.path import paths
