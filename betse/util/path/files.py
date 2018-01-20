#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level non-directory filename facilities.

See Also
----------
:mod:`betse.util.io.iofiles`
    Low-level non-directory file content facilities.
'''

# ....................{ IMPORTS                            }....................
import os, shutil
from betse.exceptions import BetseFileException
from betse.util.io.log import logs
from betse.util.type.types import type_check
from os import path as os_path

# ....................{ EXCEPTIONS                         }....................
def die_if_file(pathname: str) -> None:
    '''
    Raise an exception if the passed path is an existing non-directory file
    *after* following symbolic links.

    Raises
    ----------
    BetseFileException
        If any passed path is such a file.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    # If this file exists...
    if is_file(pathname):
        # Avoid circular import dependencies.
        from betse.util.path import paths

        # Raise a human-readable exception.
        raise BetseFileException(
            'Path "{}" already an existing {}.'.format(
                pathname, paths.get_type_label(pathname)))


def die_unless_file(*pathnames: str) -> None:
    '''
    Raise an exception unless all passed paths are existing non-directory files
    *after* following symbolic links.

    Raises
    ----------
    BetseFileException
        If any passed path is *NOT* such a file.

    See Also
    ----------
    :func:`is_file`
        Further details.
    '''

    # If at least one passed path is *NOT* an existing non-directory file...
    if not is_file(*pathnames):
        # For each such path...
        for pathname in pathnames:
            # If this path *NOT* such a file, raise an exception.
            if not is_file(*pathnames):
                raise BetseFileException(
                    'File "{}" not found or unreadable.'.format(pathname))


def join_and_die_unless_file(*pathnames: str) -> str:
    '''
    Pathname of the file produced by joining (i.e., concatenating) the passed
    pathnames with the directory separator specific to the current platform if
    this file exists *or* raise an exception otherwise (i.e., if this file does
    *not* exist).

    See Also
    -----------
    :func:`betse.util.path.pathnames.join`
    :func:`die_unless_file`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.pathnames import join

    # Filename producing by joining these pathnames.
    filename = join(*pathnames)

    # If this file does *NOT* exist, raise an exception.
    die_unless_file(filename)

    # Return this dirname.
    return filename

# ....................{ EXCEPTIONS ~ writable              }....................
def die_unless_exists_writable(pathname: str) -> None:
    '''
    Raise an exception unless the passed path is an existing writable
    non-directory file *after* following symbolic links.

    Raises
    ----------
    BetseFileException
        If this path is *NOT* such a file.

    See Also
    ----------
    :func:`is_writable`
        Further details.
    '''

    # If this file does not exist or does but is unwritable, raise an exception.
    if not is_exists_writable(pathname):
        raise BetseFileException(
            'File "{}" not found or unwritable.'.format(pathname))

# ....................{ TESTERS                            }....................
@type_check
def is_file(*pathnames: str) -> bool:
    '''
    ``True`` only if all passed paths are existing non-directory files *after*
    following symbolic links.

    If any such path does *not* exist, this function silently ignores this path
    rather than raising an exception.

    Versus :func:`path.isfile`
    ----------
    This function intrinsically differs from the standard :func:`path.isfile`
    function. While the latter returns ``True`` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for *all* non-directory files
    regardless of whether these files are special or not.

    **Why?** Because this function complies with POSIX semantics, whereas
    :func:`path.isfile` does *not*. The specialness of non-directory files is
    usually irrelevant; in general, it only matters whether these files are
    directories or not. For example, the external command ``rm`` removes only
    non-directory files (regardless of specialness) while the external command
    ``rmdir`` removes only empty directories.

    Parameters
    ----------
    pathnames: tuple[str]
        Tuple of all paths to be tested.

    Returns
    ----------
    bool
        ``True`` only if all such paths are existing non-directory files *after*
        following symbolic links.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths

    # Return true only if...
    return any(
        # This path both exists and is *NOT* a directory...
        paths.is_path(pathname) and not dirs.is_dir(pathname)
        # For each such path.
        for pathname in pathnames
    )


def is_executable(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an **executable non-directory file**
    (i.e., file with the execute bit enabled) *after* following symbolic links.

    For generality, this function does *not* raise an exception if this path
    does not exist or does but is an existing directory.
    '''

    return is_file(pathname) and os.access(pathname, os.X_OK)


def is_exists_writable(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an **existing writable non-directory
    file** (i.e., file with the write bit enabled for the current user or user
    group) *after* following symbolic links.

    For generality, this function does *not* raise an exception if this path
    does not exist or does but is an existing directory.

    Caveats
    ----------
    For disambiguity, this function is intentionally *not* named
    ``is_writable()``. Why? Because this function returns ``False`` when this
    path does *not* currently exist but the current user has sufficient
    permissions to create and hence write this file.
    '''

    return is_file(pathname) and os.access(pathname, os.W_OK)

# ....................{ TESTERS ~ symlink                  }....................
@type_check
def is_symlink(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an existing symbolic link *before*
    following symbolic links.

    This function does *not* raise an exception if this path does not exist.
    '''

    return os_path.islink(pathname)


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
    return is_symlink(pathname) and os_path.exists(pathname)

# ....................{ GETTERS                            }....................
@type_check
def get_size(filename: str) -> int:
    '''
    Size in bytes of the contents of the file with the passed filename.
    '''

    # If this file does *NOT* exist, raise an exception.
    die_unless_file(filename)

    # Return this filesize.
    return os_path.getsize(filename)

# ....................{ COPIERS                            }....................
@type_check
def copy(
    # Mandatory parameters.
    src_filename: str,
    trg_filename: str,

    # Optional parameters.
    is_overwritable: bool = False,
) -> None:
    '''
    Copy the passed source file to the passed target file or directory.

    If the source file is a symbolic link, this link (rather than its
    transitive target) will be copied and hence preserved.

    The target file will be copied in a manner maximally preserving metadata
    (e.g., owner, group, permissions, times, extended file system attributes).
    If the target file is a directory, the basename of the source file will be
    appended to this directory -- much like the standard ``cp`` POSIX command.

    Parameters
    ----------
    src_filename : str
        Absolute or relative path of the source file to be copied from.
    trg_filename : str
        Absolute or relative path of the target file to be copied to.
    is_overwritable : optional[bool]
        If this target file already exists and this boolean is ``True``, this
        file is silently overwritten; else, an exception is raised. Defaults to
        ``False``.

    Raises
    ----------
    BetseFileException
        If either the source file does not exist *or* the target file already
        exists.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, paths, pathnames

    # Log this copy.
    logs.log_debug('Copying file: %s -> %s', src_filename, trg_filename)

    # If this source file does *NOT* exist, raise an exception.
    die_unless_file(src_filename)

    # If this target file is a directory, append the basename of the passed
    # source file to this directory -- much like the "cp" POSIX command.
    if dirs.is_dir(trg_filename):
        trg_filename = pathnames.join(
            trg_filename, pathnames.get_basename(src_filename))

    # If this target file already exists but is *NOT* overwritable, raise an
    # exception.
    if not is_overwritable:
        paths.die_if_path(trg_filename)

    # Create the directory containing this target file *BEFORE* calling the
    # shutil.copy2() function, which assumes this directory to exist.
    dirs.make_parent_unless_dir(trg_filename)

    # Perform this copy in a manner preserving metadata and symbolic links.
    shutil.copy2(src_filename, trg_filename, follow_symlinks=False)

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
