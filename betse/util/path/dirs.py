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
import os, shutil

# ....................{ EXCEPTIONS                         }....................
def die_unless_dir(*dirnames) -> None:
    '''
    Raise an exception unless all passed directories exist.
    '''
    for dirname in dirnames:
        if not is_dir(dirname):
            raise BetseExceptionDir(
                'Directory "{}" not found or not a readable directory.'.format(
                    dirname))

def die_unless_parent_dir(pathname: str) -> None:
    '''
    Raise an exception unless the parent directory of the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    die_unless_dir(paths.get_dirname(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_dir(*dirnames) -> None:
    '''
    Raise an exception if any of the passed directories exist.
    '''
    for dirname in dirnames:
        if is_dir(dirname):
            raise BetseExceptionDir(
                'Directory "{}" already exists.'.format(dirname))

# ....................{ TESTERS                            }....................
def is_dir(dirname: str) -> bool:
    '''
    True if the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    return path.isdir(dirname)

# ....................{ LISTERS                            }....................
def list_basenames(dirname: str) -> list:
    '''
    Get a list of the basenames of all files and subdirectories in the passed
    directory.
    '''
    die_unless_dir(dirname)
    return os.listdir(dirname)

# ....................{ MAKERS                             }....................
def make_unless_dir(dirname: str) -> None:
    '''
    Create the passed directory if such directory does *not* already exist.

    All nonexistent parents of such directory will also be recursively created,
    mimicking the action of the standard `mkdir -p` shell command.
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

def make_parent_unless_dir(*pathnames) -> None:
    '''
    Create the parent directory of each passed path for any such directory that
    does *not* already exist.

    All nonexistent parents of each such directory will also be recursively
    created, mimicking the action of the standard `mkdir -p` shell command.
    '''
    # Canonicalize each pathname *BEFORE* attempting to get its dirname.
    # Relative pathnames do *NOT* have sane dirnames (e.g., the dirname for a
    # relative pathname "metatron" is the empty string) and hence *MUST* be
    # converted to absolute pathnames first.
    for pathname in pathnames:
        make_unless_dir(paths.get_dirname(paths.canonicalize(pathname)))

# ....................{ COPIERS                            }....................
def copy_into_target_dir(dirname_source: str, dirname_target: str) -> None:
    '''
    Recursively copy the passed source directory to a subdirectory of the passed
    target directory having the same basename as such source directory.

    See Also
    ----------
    copy()
        For further details.

    Examples
    ----------
    >>> from betse.util.path import dirs
    >>> dirs.copy_into_target_dir('/usr/src/linux/', '/tmp/')
    >>> dirs.is_dir('/tmp/linux/')
    True
    '''
    assert isinstance(dirname_target, str),\
        '"{}" not a string.'.format(dirname_target)
    assert len(dirname_target), 'Target dirname empty.'

    # Perform such copy.
    basename_source = paths.get_basename(dirname_source)
    copy(dirname_source, paths.join(dirname_target, basename_source))

def copy(dirname_source: str, dirname_target: str) -> None:
    '''
    Recursively copy the passed source to target directory.

    All nonexistent parents of the target directory will be recursively created,
    mimicking the action of the `mkdir -p` shell command. All symbolic links in
    the source directory will be preserved (i.e., copied as is rather than their
    transitive targets copied in their stead).

    If the source directory does not exist *or* the target directory already
    exists, an exception will be raised.
    '''
    assert isinstance(dirname_source, str),\
        '"{}" not a string.'.format(dirname_source)
    assert isinstance(dirname_target, str),\
        '"{}" not a string.'.format(dirname_target)
    assert len(dirname_source), 'Source dirname empty.'
    assert len(dirname_target), 'Target dirname empty.'

    # Log such copy.
    loggers.log_info(
        'Copying directory "%s" to "%s".', dirname_source, dirname_target)

    # Raise an exception if the target directory already exists. While we could
    # defer to the exception raised by the shutil.copytree() function for such
    # case, such exception's message erroneously refers to such directory as a
    # file and is hence best avoided: e.g.,
    #
    #     [Errno 17] File exists: 'sample_sim'
    die_if_dir(dirname_target)

    # Perform such copy.
    shutil.copytree(
        src = dirname_source,
        dst = dirname_target,
        symlinks = True,
    )

# --------------------( WASTELANDS                         )--------------------
# def copy_unless_dir(dirname_source: str, dirname_target: str) -> None:
#     '''
#     Recursively copy the passed source to target directory unless the latter
#     already exists, in which case this function reduces to a noop.
#
#     See Also
#     ----------
#     copy()
#         For further details.
#     '''
#     if not is_dir(dirname_target):
#         copy(dirname_source, dirname_target)

#FUXME: Replace all existing calls to os.makedirs() by calls to such functions.
# from betse.util.path import paths
