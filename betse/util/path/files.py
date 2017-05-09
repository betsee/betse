#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level non-directory filename facilities.

See Also
----------
:mod:`betse.util.io.files`
    Low-level non-directory file content facilities.
'''

# ....................{ IMPORTS                            }....................
import os, shutil
from os import path
from betse.exceptions import BetseFileException
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_file(pathname: str) -> None:
    '''
    Raise an exception unless the passed path is a non-directory file *after*
    following symbolic links.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    if not is_file(pathname):
        raise BetseFileException(
            'File "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed path is a non-directory file *after*
    following symbolic links.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    if is_file(pathname):
        raise BetseFileException('File "{}" already exists.'.format(pathname))


def die_if_special(pathname: str) -> None:
    '''
    Raise an exception if the passed path is a special file.

    See Also
    ----------
    :func:`is_special`
        Further details.
    '''

    if is_special(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import paths

        raise BetseFileException(
            'File "{}" already an existing {}.'.format(
                pathname, paths.get_type_label(pathname)))

# ....................{ TESTERS                            }....................
def is_file(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is a non-directory file _after_ following
    symbolic links.

    This function does *not* raise an exception if this path does not exist.

    Versus `path.isfile()`
    ----------
    This function intrinsically differs from the standard :func:`path.isfile`
    function. While the latter returns ``True`` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for _all_ non-directory files
    regardless of whether these files are special or not.

    **Why?** Because this function complies with POSIX semantics, whereas
    :func:`path.isfile` does *not*. The specialness of non-directory files is
    usually irrelevant; in general, it only matters whether these files are
    directories or not. For example, the external command `rm` removes only
    non-directory files (regardless of specialness) while the external command
    `rmdir` removes only empty directories.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # This path is a file if this path both exists and is *NOT* a directory.
    return paths.is_path(pathname) and not dirs.is_dir(pathname)


def is_file_executable(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an **executable non-directory file**
    (i.e., file with the execute bit enabled) _after_ following symbolic links.

    This function does *not* raise an exception if this path does not exist.
    '''

    # This path is an executable file if this path is an existing file with the
    # executable bit enabled.
    return is_file(pathname) and os.access(pathname, os.X_OK)


#FIXME: Given that this function resides in the "files" submodule, shouldn't
#this function return False rather than True when passed a directory pathname?
#On the other hand, if this function is currently working as intended, this
#function should be shifted into the "betse.util.path.paths" submodule.
def is_special(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an existing **special file** (e.g.,
    directory, device node, socket, symbolic link).

    This function does *not* raise an exception if this path does not exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # True if this path exists and...
    return paths.is_path(pathname) and (
        # ...is either a symbolic link *OR* neither a regular file nor symbolic
        # link to such a file. In the latter case, predicate logic guarantees
        # this file to *NOT* be a symbolic link, thus reducing this test to:
        # "...is either a symbolic link *OR* not a regular file."
        is_symlink(pathname) or not path.isfile(pathname))

# ....................{ TESTERS ~ symlink                  }....................
@type_check
def is_symlink(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an existing symbolic link *before*
    following symbolic links.

    This function does *not* raise an exception if this path does not exist.
    '''

    return path.islink(pathname)


@type_check
def is_symlink_valid(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an existing **non-dangling symbolic link**
    (i.e., symbolic link whose target also exists) *before* following symbolic
    links.

    This function does *not* raise an exception if this path does not exist.
    '''

    # Call path.exists() rather than path.lexists(), as the latter returns True
    # for dangling symbolic links.
    #
    # This is why human-readable function names is a good thing, people.
    return is_symlink(pathname) and path.exists(pathname)

# ....................{ COPIERS                            }....................
@type_check
def copy(filename_source: str, filename_target: str) -> None:
    '''
    Copy the passed source file to the passed target file or directory.

    If the source file is a symbolic link, this link (rather than its
    transitive target) will be copied and hence preserved.

    The target file will be copied in a manner maximally preserving metadata
    (e.g., owner, group, permissions, times, extended file system attributes).
    If the target file is a directory, the basename of the source file will be
    appended to this directory -- much like the standard `cp` POSIX command.

    If either the source file does not exist *or* the target file already
    exists, an exception will be raised.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Log this copy.
    logs.log_debug(
        'Copying file "%s" to "%s".', filename_source, filename_target)

    # Raise an exception unless the source file exists.
    die_unless_file(filename_source)

    # If the target file is a directory, append the basename of the passed
    # source file to this directory -- much like the "cp" POSIX command.
    if dirs.is_dir(filename_target):
        filename_target = paths.join(
            filename_target, paths.get_basename(filename_source))

    # Raise an exception if the target file already exists.
    paths.die_if_path(filename_target)

    # Perform this copy in a manner preserving metadata and symbolic links.
    shutil.copy2(filename_source, filename_target, follow_symlinks=False)

# ....................{ REMOVERS                           }....................
@type_check
def remove_if_found(filename: str) -> None:
    '''
    Remove the passed non-directory file if this file currently exists.

    If this file does *not* currently exist, this function reduces to a noop.
    For safety, this function removes this file atomically; in particular, this
    file's existence is *not* explicitly tested for.
    '''

    # Log this removal if the subsequent removal attempt is likely to actually
    # remove a file. Due to race conditions with other threads and processes,
    # this file could be removed after this test succeeds but before the removal
    # is performed. Since this is largely ignorable, the worst case is an
    # extraneous log message.
    if is_file(filename):
        logs.log_debug('Removing file: %s', filename)

    # Remove this file atomically. To avoid race conditions with other threads
    # and processes, this operation is *NOT* embedded in an explicit test for
    # file existence. Instead, the Pythonic Way is embraced.
    try:
        os.remove(filename)
    # If this file does *NOT* exist, ignore this exception.
    except FileNotFoundError:
        pass


@type_check
def remove(filename: str) -> None:
    '''
    Remove the passed non-directory file.
    '''

    # Log this removal.
    logs.log_debug('Removing file: %s', filename)

    # Raise an exception unless this file exists.
    die_unless_file(filename)

    # Remove this file. Note that the os.remove() and os.unlink() functions are
    # identical. (That was a tad silly, Guido.)
    os.remove(filename)
