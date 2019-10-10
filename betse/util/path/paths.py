#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level path facilities general to all path types -- namely, both directories
and non-directory files.
'''

# ....................{ IMPORTS                           }....................
import os, shutil
from betse.exceptions import BetsePathException
from betse.util.io.log import logs
from betse.util.type.types import type_check, IterableTypes, NumericSimpleTypes
from os import path as os_path

# ....................{ EXCEPTIONS                        }....................
def die_if_path(*pathnames: str) -> None:
    '''
    Raise an exception if any of the paths with the passed pathnames exist.

    Parameters
    ----------
    pathnames : tuple[pathnames]
        Tuple of all pathnames to be tested.

    Raises
    ----------
    BetsePathException
        If any such path exists.
    '''

    for pathname in pathnames:
        if is_path(pathname):
            raise BetsePathException(
                'Path "{}" already exists.'.format(pathname))


def die_unless_path(*pathnames: str) -> None:
    '''
    Raise an exception unless all paths with the passed pathnames exist.

    Parameters
    ----------
    pathnames : tuple[pathnames]
        Tuple of all pathnames to be tested.

    Raises
    ----------
    BetsePathException
        If any such path does *not* exist.
    '''

    for pathname in pathnames:
        if not is_path(pathname):
            raise BetsePathException(
                'Path "{}" not found or unreadable.'.format(pathname))


@type_check
def path_or_die(pathname: str) -> str:
    '''
    Passed pathname if a path with this pathname exists *or* raise an exception
    otherwise (i.e., if this path does *not* exist).

    Parameters
    ----------
    pathname : str
        Absolute or relative pathname of the path to be validated.

    Returns
    ----------
    str
        Passed pathname if this path exists.

    Raises
    ----------
    BetsePathException
        If this path does *not* exist.
    '''

    # If this path is not found, fail.
    die_unless_path(pathname)

    # Else, this path exists. Return this path's pathname.
    return pathname

# ....................{ EXCEPTIONS ~ type                 }....................
def die_if_special(pathname: str) -> None:
    '''
    Raise an exception if the passed path is an existing special file.

    Raises
    ----------
    BetsePathException
        If this path is *not* an existing special file.

    See Also
    ----------
    :func:`is_special`
        Further details.
    '''

    # If this special file exists, raise a human-readable exception.
    if is_special(pathname):
        raise BetsePathException(
            'Path "{}" already an existing {}.'.format(
                pathname, get_type_label(pathname)))

# ....................{ EXCEPTIONS ~ permissions          }....................
def die_unless_writable(pathname: str) -> None:
    '''
    Raise an exception unless the passed path is writable by the current user.

    Raises
    ----------
    BetsePathException
        If this path is *not* writable by this user.

    See Also
    ----------
    :func:`is_special`
        Further details.
    '''

    # If this path is unwritable, raise a human-readable exception.
    if not is_writable(pathname):
        raise BetsePathException(
            'Path "{}" unwritable by current user.'.format(pathname))

# ....................{ TESTERS                           }....................
@type_check
def is_path(pathname: str) -> bool:
    '''
    ``True`` only if the passed path exists *without* following this path if
    this path is an existing symbolic link.

    If this path is an existing **dangling symbolic link** (i.e., symbolic link
    whose target no longer exists, also commonly referred to as a broken
    symbolic link), this function still returns ``True``.
    '''

    # Call path.lexists() rather than path.exists(), as the latter returns
    # False for dangling symbolic links -- which is entirely irrelevant for
    # most contexts. Under POSIX semantics, dangling symbolic links are
    # essential to common usage patterns and should *NOT* be discriminated
    # against here.
    return os_path.lexists(pathname)


def is_special(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is an existing **special file** (e.g.,
    directory, device node, socket, symbolic link).

    This function does *not* raise an exception if this path does not exist.
    '''

    # Avoid circular import dependencies.
    from betse.util.file import files

    # True if this path exists and...
    return is_path(pathname) and (
        # ...is either a symbolic link *OR* neither a regular file nor symbolic
        # link to such a file. In the latter case, predicate logic guarantees
        # this file to *NOT* be a symbolic link, thus reducing this test to:
        # "...is either a symbolic link *OR* not a regular file."
        files.is_symlink(pathname) or not os_path.isfile(pathname))

# ....................{ TESTERS ~ permissions             }....................
@type_check
def is_readable(pathname: str) -> bool:
    '''
    ``True`` only if the passed path both exists *and* is readable by the
    current user.
    '''

    return is_path(pathname) and os.access(pathname, os.R_OK)


@type_check
def is_writable(pathname: str) -> bool:
    '''
    ``True`` only if the passed pathname is **writable** (i.e., either an
    existing path with the write bit enabled for the current user or user group
    *or* a non-existing path whose parent directory has that bit enabled)
    *after* following symbolic links.

    For generality, this function does *not* raise an exception if this path
    does not exist or does but is an existing directory.

    See Also
    ----------
    :func:`files.is_exists_writable`
        Similar tester adding a file-specific existence requirement.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, pathnames

    # If this path exists, return true only if this path is writable.
    if is_path(pathname):
        return os.access(pathname, os.W_OK)
    # Else, this path does *NOT* exist. In this case...
    else:
        # Absolute or relative dirname of this path's parent directory.
        dirname = pathnames.get_dirname(pathname)

        # Return true only if this directory both exists and is writable.
        return dirs.is_dir(dirname) and os.access(dirname, os.W_OK)

# ....................{ TESTERS ~ mtime : recursive        }....................
@type_check
def is_mtime_recursive_older_than_paths(
    pathname: str, pathnames: IterableTypes) -> bool:
    '''
    ``True`` only if the **recursive mtime** (i.e., recursively calculated
    modification time) of the passed path is strictly less than that of at
    least one path in the passed iterable.

    For efficiency, this function short-circuits at the first path in this
    iterable newer than than this path. Callers are thus recommended to pass an
    iterable pre-sorted such that:

    * Non-directory files (whose recursive mtime reduces to an efficient
      constant-time operation) are ordered *before* directories (whose recursive
      mtime requires performing an inefficient filesystem walk).
    * Directories expected to contain fewer files and subdirectories are ordered
      *before* directories expected to contain more files and subdirectories.

    This function is functionally equivalent to (but substantially faster in
    both the worst and common case than) the following conditional logic:

    .. code:: python

       get_mtime_recursive(pathname) < get_mtime_recursive_newest(pathnames)

    Parameters
    ----------
    pathname : str
        Absolute or relative pathname of the path to be tested.
    pathnames: IterableTypes[str]
        Iterable of the absolute or relative pathnames of paths to test this
        path against.

    Returns
    ----------
    bool
        ``True`` only if this path's recursive mtime is strictly less than that
        of at least one path in this iterable.
    '''

    # Recursive mtime of this path.
    path_mtime_recursive = get_mtime_recursive(pathname)

    # Return True only if...
    return any(
        # This mtime is less than that of at least one path in this iterable.
        path_mtime_recursive < get_mtime_recursive(pathname)
        for pathname in pathnames
    )

# ....................{ GETTERS                            }....................
@type_check
def get_type_label(pathname: str) -> str:
    '''
    Lowercase human-readable label describing the type of the passed existing
    path.

    For brevity, this label comprises one lowercase noun optionally preceded by
    a lowercase adjective unambiguously identifying this path's type:

    * ``symbolic link`` for an existing symbolic link, regardless of whether
      that link is **dangling** (i.e., refers to a non-existent path) or not.
    * ``directory`` for an existing directory.
    * ``special file`` for all other **special files** (e.g., device node, named
      pipe, socket).
    * ``file`` for an existing regular (i.e., non-special) file).

    If this path does *not* exist, an exception is raised.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, files

    # Raise an exception unless this path exists.
    die_unless_path(pathname)

    # Get this type.
    if files.is_file(pathname):
        if files.is_symlink(pathname):
            return 'symbolic link'
        elif is_special(pathname):
            return 'special file'
        else:
            return 'file'
    elif dirs.is_dir(pathname):
        return 'directory'

    # Else, this existing path is neither a file or directory and is hence a
    # shambling thing that should not be. Panic. Space-time continuity itself is
    # now in question.
    raise BetsePathException(
        'Path "{}" type unrecognized.'.format(pathname))

# ....................{ GETTERS ~ mtime : non-recursive    }....................
@type_check
def get_mtime_nonrecursive(pathname: str) -> NumericSimpleTypes:
    '''
    **Non-recursive mtime** (i.e., non-recursively retrieved modification time) in
    in seconds of the passed path.

    If this path is:

    * A directory, this is the most recent time at which the non-recursive
      contents of this directory were created, modified, or deleted. Since this
      function operates non-recursively, this considers *only* this directory
      itself and all immediate subdirectories and files in this directory.
    * A file, this is the most recent time at which the contents of this file
      were created or modified.
    * A symbolic link, this is the most recent time at which this symbolic link
      was created. Hence, this function does *not* follow symbolic links.

    This time does *not* consider the most recent time at which the metadata
    (e.g., permissions) of this path was created or modified, which the
    :func:`get_ctime` function provides.

    Parameters
    -----------
    pathname : str
        Relative or absolute path to inspect.

    Returns
    -----------
    NumericSimpleTypes
        Mtime in seconds of this path as either an integer or float, depending
        on the boolean returned by the platform-specific
        :func:`os.stat_float_times` function.
    '''

    return os_path.getmtime(pathname)

# ....................{ GETTERS ~ mtime : recursive        }....................
@type_check
def get_mtime_recursive(pathname: str) -> NumericSimpleTypes:
    '''
    **Recursive mtime** (i.e., recursively calculated modification time) in
    seconds of the passed path.

    If this path is:

    * A directory, this is the most recent time at which the recursive contents
      of this directory were created, modified, or deleted. Since this function
      operates recursively, this considers this directory itself and all
      transitive subdirectories and files of this directory.
    * A file, this is the most recent time at which the contents of this file
      were created or modified.
    * A symbolic link, this is the most recent time at which this symbolic link
      was created. Hence, this function does *not* follow symbolic links.

    Parameters
    -----------
    pathnames : IterableTypes[str]
        Iterable of all relative and absolute paths to inspect.

    Returns
    -----------
    NumericSimpleTypes
        Mtime in seconds of the most recent path as either an integer or float,
        depending on the boolean returned by the platform-specific
        :func:`os.stat_float_times` function.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Return, if this path is a:
    return (
        # * Directory, the mtime of the most recent path in this directory.
        dirs.get_mtime_recursive_newest(pathname) if dirs.is_dir(pathname) else
        # * Non-directory file, the mtime of this file as is.
        get_mtime_nonrecursive(pathname)
    )


@type_check
def get_mtime_recursive_newest(pathnames: IterableTypes) -> NumericSimpleTypes:
    '''
    **Recursive mtime** (i.e., recursively calculated modification time) in
    seconds of the most recent passed path.

    Parameters
    -----------
    pathnames : IterableTypes[str]
        Iterable of all relative and absolute paths to be inspected.

    Returns
    -----------
    NumericSimpleTypes
        Mtime in seconds of the most recent path as either an integer or float,
        depending on the boolean returned by the platform-specific
        :func:`os.stat_float_times` function.

    See Also
    -----------
    :func:`get_mtime_recursive`
        Further details.
    '''

    # Return the maximum of a generator expression yielding the recursive mtimes
    # of all passed paths.
    return max(get_mtime_recursive(pathname) for pathname in pathnames)

# ....................{ MOVERS                             }....................
@type_check
def move_path(src_pathname: str, trg_pathname: str) -> None:
    '''
    Move the source path with the passed pathname to the target path with the
    passed pathname.

    This function moves this source path in a manner that maximally preserves
    metadata (e.g., owner, group, permissions, times, extended file system
    attributes). Likewise, if this source path is a symbolic link, this
    function moves (and hence preserves) this link rather than its transitive
    target.

    Parameters
    ----------
    src_pathname: str
        Absolute or relative pathname of the source path to be copied from.
    trg_pathname: str
        Absolute or relative pathname of the target path to be copied to.

    Raises
    ----------
    BetsePathException
        If either:

        * The source path does *not* exist.
        * The target path already exists.
    '''

    # Log this move.
    logs.log_debug(
        'Moving %s: %s -> %s',
        get_type_label(src_pathname), src_pathname, trg_pathname)

    # If this source path does *NOT* exist, raise an exception.
    die_unless_path(src_pathname)

    # If the target path already exists, raise an exception.
    #
    # Note that this is essential to sane cross-platform semantics. The
    # shutil.move() function is implemented in terms of os.rename(), which
    # overwrites the passed path if that path is a non-directory file
    # contextually depending on platform. (In most cases, that's bad.)
    die_if_path(trg_pathname)

    # Move this source to target path.
    shutil.move(src_pathname, trg_pathname)

# ....................{ REMOVERS                          }....................
@type_check
def remove_path(pathname: str) -> None:
    '''
    Recursively remove the passed path in a reasonable manner (e.g., *not*
    following symbolic links outside this path if this path is a directory).

    Caveats
    ----------
    Recursive directory removal is an inherently dangerous operation. If this
    path is a directory rather than a file, this function (in order):

    1. Notifies the end user with a logged warning.
    1. Waits several seconds, enabling sufficiently aware end users to jam the
       panic button.
    1. Recursively removes this directory.

    Parameters
    ----------
    pathname: str
        Absolute or relative pathname of the source path to be copied from.

    Raises
    ----------
    BetsePathException
        If this path does *not* exist.

    See Also
    ----------
    :func:`betse.util.path.dirs.remove_dir`
    :func:`betse.util.path.files.remove_file`
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs, files

    # If this path does *NOT* exist, raise an exception.
    die_unless_path(pathname)

    # If this path is a directory, remove this directory.
    if dirs.is_dir(pathname):
        dirs.remove_dir(pathname)
    # Else, remove this file.
    else:
        files.remove_file(pathname)
