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
import os, shutil
from betse.exceptions import BetseExceptionDir
from betse.util.io import loggers
from betse.util.type import types
from os import path

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_dir(*dirnames) -> None:
    '''
    Raise an exception unless all passed directories exist.
    '''
    for dirname in dirnames:
        if not is_dir(dirname):
            raise BetseExceptionDir(
                'Directory "{}" not found or unreadable.'.format(dirname))


def die_unless_parent_dir(pathname: str) -> None:
    '''
    Raise an exception unless the parent directory of the passed path exists.
    '''
    # Avoid circular import dependencies.
    from betse.util.path import paths
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
    `True` if the passed directory exists.
    '''
    assert types.is_str_nonempty(dirname),\
        types.assert_not_str_nonempty(dirname, 'Dirname')
    return path.isdir(dirname)

# ....................{ GETTERS                            }....................
def get_current_dirname() -> str:
    '''
    Get the **current working dirname** (i.e., absolute path of the current
    working directory (CWD)).

    Unless subsequently changed, this is the absolute path of the directory from
    which `betse` was first run.
    '''
    return os.getcwd()

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
    Create the passed directory if this directory does *not* already exist.

    All nonexistent parents of this directory will also be recursively created,
    mimicking the action of the standard `mkdir -p` shell command.
    '''
    assert types.is_str_nonempty(dirname), (
        types.assert_not_str_nonempty(dirname, 'Dirname'))

    # If this directory does *NOT* already exist, create this directory. To
    # support logging, this condition is explicitly tested for. To avoid race
    # conditions (e.g., in the event this directory is created between testing
    # and creating this directory), we preserve the makedirs() keyword argument
    # "exist_ok = True" below.
    if not is_dir(dirname):
        # Log this creation.
        loggers.log_debug('Creating directory "%s".', dirname)

        # Create this directory if still needed.
        os.makedirs(dirname, exist_ok = True)


def canonicalize_and_make_unless_dir(dirname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed
    directory _and_ create this directory if this directory does *not* already
    exist.

    This convenience function simply passes this directory to the
    `paths.canonicalize()` and `dirs.make_unless_dir()` functions (in that
    order) and then returns this directory.

    See Also
    -----------
    make_unless_dir()
        Further details.
    '''
    # Avoid circular import dependencies.
    from betse.util.path import paths
    dirname = paths.canonicalize(dirname)
    make_unless_dir(dirname)
    return dirname


def make_parent_unless_dir(*pathnames) -> None:
    '''
    Create the parent directory of each passed path for any such directory that
    does *not* already exist.

    All nonexistent parents of each such directory will also be recursively
    created, mimicking the action of the standard `mkdir -p` shell command.
    '''
    # Avoid circular import dependencies.
    from betse.util.path import paths

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
    assert types.is_str_nonempty(dirname_source),\
        types.assert_not_str_nonempty(dirname_source, 'Source dirname')
    assert types.is_str_nonempty(dirname_target),\
        types.assert_not_str_nonempty(dirname_target, 'Target dirname')

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Copy us up the directory bomb.
    basename_source = paths.get_basename(dirname_source)
    copy(dirname_source, paths.join(dirname_target, basename_source))


def copy(dirname_source: str, dirname_target: str) -> None:
    '''
    Recursively copy the passed source to target directory.

    All nonexistent parents of the target directory will be recursively created,
    mimicking the action of the `mkdir -p` shell command. All symbolic links in
    the source directory will be preserved (i.e., copied as is rather than their
    transitive targets copied instead).

    If either the source directory does not exist *or* the target directory
    already exists, an exception will be raised.
    '''
    assert types.is_str_nonempty(dirname_source),\
        types.assert_not_str_nonempty(dirname_source, 'Source dirname')
    assert types.is_str_nonempty(dirname_target),\
        types.assert_not_str_nonempty(dirname_target, 'Target dirname')

    # Log such copy.
    loggers.log_info(
        'Copying directory "%s" to "%s".', dirname_source, dirname_target)

    # Raise an exception unless the source directory exists.
    die_unless_dir(dirname_source)

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
