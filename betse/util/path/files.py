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
from betse.exceptions import BetseExceptionFile
from betse.util.io import loggers
from os import path
import os, shutil

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

# ....................{ COPIERS                            }....................
def copy(filename_source: str, filename_target: str) -> None:
    '''
    Copy the passed source to target non-directory file.

    Such file will be copied in a manner maximally preserving metadata (e.g.,
    owner, group, permissions, times, extended file system attributes). If such
    source file is a symbolic link, such link (rather than the transitive
    target of such link) will be copied.
    '''
    assert isinstance(filename_source, str),\
        '"{}" not a string.'.format(filename_source)
    assert isinstance(filename_target, str),\
        '"{}" not a string.'.format(filename_target)

    # Fail unless such source file exists.
    die_unless_found(filename_source)

    # Log such copy.
    loggers.log_info(
        'Copying file "%s" to "%s".', filename_source, filename_target)

    # Copy such source file, preserving metadata and symbolic links.
    shutil.copy2(filename_source, filename_target, follow_symlinks = False)

# ....................{ REMOVERS                           }....................
def remove(filename: str) -> None:
    '''
    Remove the passed non-directory file.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)

    # Fail unless such file exists.
    die_unless_found(filename)

    # Log such removal.
    loggers.log_info('Removing file "%s".', filename)

    # Remove such file. Note that the os.remove() and os.unlink() functions are
    # identical. (That was silly, Guido.)
    os.remove(filename)

# --------------------( WASTELANDS                         )--------------------
    # Such file will be copied in a manner preserving some but *not* all metadata,
    # in accordance with standard POSIX behaviour. Specifically, the permissions
    # but *not* owner, group, or times of such file

    # Such file will be copied in a manner maximally preserving metadata (e.g.,
    # owner, group, permissions, times, extended file system attributes
    # If such source file is a symbolic link, such link rather than the target
    # file of such link will be copied.
    # If such file does *not* exist, an exception is raised.
