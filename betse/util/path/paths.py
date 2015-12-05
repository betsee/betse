#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level path facilities specific to neither directories nor non-directory
files.

This module is named `paths` rather than `path` to avoid conflict with the stock
`path` module of the `os` package.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionPath
from betse.util.io import loggers
from betse.util.type import types
from os import path
import shutil

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_path(*pathnames) -> None:
    '''
    Raise an exception unless all passed paths exist.
    '''
    for pathname in pathnames:
        if not is_path(pathname):
            raise BetseExceptionPath(
                'Path "{}" not found or unreadable.'.format(pathname))


def die_if_dirname_empty(pathname: str) -> None:
    '''
    Raise an exception if the passed path has no dirname and hence contains no
    directory separators.

    See Also
    ----------
    `is_dirname_empty()`
        For further details.
    '''
    if is_dirname_empty(pathname):
        raise BetseExceptionPath(
            'Path "{}" contains no directory separators.'.format(pathname))

# ....................{ EXCEPTIONS ~ if                    }....................
def die_if_path(*pathnames) -> None:
    '''
    Raise an exception if any passed path exists.
    '''
    for pathname in pathnames:
        if is_path(pathname):
            raise BetseExceptionPath(
                'Path "{}" already exists.'.format(pathname))

# ....................{ TESTERS                            }....................
def is_path(pathname: str) -> bool:
    '''
    `True` if the passed path exists.

    If such path is an existing **broken symbolic link** (i.e., a symbolic link
    whose target no longer exists), this function still returns True.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')

    # Call path.lexists() rather than path.exists(), as the latter returns False
    # for dangling symbolic links -- which is entirely irrelevant to most
    # contexts. Under POSIX semantics, dangling symbolic links are essential to
    # common usage patterns and should *NOT* be discriminated against here.
    return path.lexists(pathname)

def is_pathname(pathname: str) -> bool:
    '''
    `True` if the passed pathname is a valid pathname for the filesystem on
    which that pathname either does or would hypothetically reside; `False`
    otherwise.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')

    # The following call only raises exceptions when passed invalid pathnames.
    try:
        path.exists(pathname)
    except Exception:
        return False
    return True

# ....................{ TESTERS                            }....................
def is_absolute(pathname: str) -> bool:
    '''
    `True` if the passed path is absolute.

    The definition of "absolute" depends on the current operating system. Under:

    * POSIX-compatible systems, absolute paths are prefixed by `/`.
    * Microsoft Windows, absolute paths are prefixed by an optional drive
      indicator (e.g., `C:`) followed by `\`.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.isabs(pathname)


def is_relative(pathname: str) -> bool:
    '''
    `True` if the passed path is relative.

    The definition of "relative" depends on the current operating system. Under:

    * POSIX-compatible systems, relative paths are _not_ prefixed by `/`.
    * Microsoft Windows, relative paths are _not_ prefixed by an optional drive
      indicator (e.g., `C:`) followed by `\`.
    '''
    return not is_absolute(pathname)

# ....................{ TESTERS ~ pathname                 }....................
def is_dirname_empty(pathname: str) -> bool:
    '''
    `True` if the passed pathname is a **pure basename** (i.e., contains no
    directory separators and hence no directory components).
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.sep not in pathname


def is_filetype(pathname: str, filetype: str) -> bool:
    '''
    `True` if the passed pathname has the passed filetype.

    Such filetype may contain arbitrarily many `.` characters, including an
    optional prefixing `.`. Regardless, this function behaves as expected.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    assert types.is_str_nonempty(filetype),\
        types.assert_nonstr_nonempty(filetype, 'Filetype')

    # Avoid circular import dependencies.
    from betse.util.type import strs

    # Test such filetype, prefixed by "." unless already prefixed.
    return strs.is_suffix(
        pathname, strs.add_prefix_unless_found(filetype, '.'))

# ....................{ GETTERS                            }....................
def get_basename(pathname: str) -> str:
    '''
    Get the **basename** (i.e., last component) of the passed path.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.basename(pathname)


def get_type_label(pathname: str) -> str:
    '''
    Get a lowercase human-readable label describing the type of the passed path.

    For brevity, such label comprises one lowercase noun optionally preceded by
    a lowercase adjective unambiguously identifying such path's type.

    If such path does *not* exist, an exception will be raised.
    '''
    # Avoid circular import dependencies.
    from betse.util.paths import dirs, files

    # Raise an exception unless such path exists.
    die_unless_path(pathname)

    # Get such type.
    if files.is_file(pathname):
        if files.is_symlink(pathname):
            return 'symbolic link'
        elif files.is_special(pathname):
            return 'special file'
        else:
            return 'file'
    elif dirs.is_dir(pathname):
        return 'directory'
    else:
        return 'path'

    # This should never happen. That means it will in several seconds.
    return None

# ....................{ GETTERS ~ dirname                  }....................
def get_dirname(pathname: str) -> str:
    '''
    Get the **dirname** (i.e., parent directory) of the passed path if such path
    has a dirname or raise an exception otherwise.
    '''
    # Raise an exception unless this path has a dirname.
    die_if_dirname_empty(pathname)

    # Get this dirname.
    dirname = get_dirname_or_empty(pathname)

    # Assert this dirname's non-emptiness. Technically, the above call *SHOULD*
    # have ensured this. This is a sufficiently critical function, however, to
    # warrant asserting this fact.
    assert len(dirname), 'Pathname "{}" dirname empty.'.format(pathname)
    return dirname


def get_dirname_or_current_dirname(pathname: str) -> str:
    '''
    Get the **dirname** (i.e., parent directory) of the passed path if such path
    has a dirname or the current working directory otherwise.
    '''
    # Avoid circular import dependencies.
    from betse.util.path import dirs
    dirname = get_dirname_or_empty(pathname)
    return dirname if dirname else dirs.get_current_dirname()


def get_dirname_or_empty(pathname: str) -> str:
    '''
    Get the **dirname** (i.e., parent directory) of the passed path if such path
    has a dirname or the empty string otherwise.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.dirname(pathname)

# ....................{ GETTERS ~ filetype                 }....................
def get_filetype(pathname: str) -> str:
    '''
    Get the **last filetype** (i.e., last `.`-prefixed substring of the
    basename *not* including that `.`) of the passed path if this path has a
    filetype or `None` otherwise.

    If this path has multiple filetypes (e.g., `odium.reigns.tar.gz`), only the
    last such filetype will be returned.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')

    # Filetype. Yes, splitext() is exceedingly poorly named.
    filetype = path.splitext(pathname)[1]

    # Strip the prefixing "." from the string returned by the prior call if this
    # path has a filetype or return None otherwise.
    return filetype[1:] if filetype else None


def get_pathname_sans_filetype(pathname: str) -> str:
    '''
    Get the passed path without last filetype (including prefixing `.`) if such
    path has a filetype *or* as is otherwise.
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.splitext(pathname)[0]

# ....................{ JOINERS                            }....................
#FIXME: According to the Python documentation, os.path.join() implicitly
#performs the following hideous operation:
#
#    If a component is an absolute path, all previous components are thrown away
#    and joining continues from the absolute path component.
#
#This is inherently dumb and we wish it to stop. Fortunately, we can! Maybe?
#It's trivial to strip leading directory separators from all passed paths except
#the first as follows:
#
#    pathnames_munged = [pathnames[0]]
#    pathnames_munged.extend(
#        pathname[1:] if path.isabs(pathname) else pathname
#        for pathname in pathnames[1:]
#    )
#
#Unfortunately, that fails to take into account the drive letter prefixing
#absolute Windows pathnames. There appears to exist a function path.splitdrive()
#doing so, but such function inefficiently returns a tuple. "Who cares about
#efficiency under Windows?" is my retort! *shrug*

def join(*pathnames) -> str:
    '''
    Join the passed pathnames on the directory separator specific to the current
    operating system.

    This is a convenience function wrapping the standard `os.path.join()`
    function _without_ adding functionality to such function -- principally to
    unify and hence simplify `import` statements in other modules.
    '''

    return path.join(*pathnames)

# ....................{ CANONICALIZERS                     }....................
def canonicalize(pathname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed path.

    Specifically (in order):

    . Perform **tilde expansion,** replacing a `~` character prefixing such path
      by the absolute path of the current user's home directory.
    . Perform **path normalization,** thus (in no particular order):
      * Collapsing redundant separators (e.g., converting `//` to `/`).
      * Converting explicit relative to absolute path components (e.g.,
        converting `../` to the name of the parent directory of such component).
      * Converting implicit relative basenames to absolute paths (e.g.,
        converting `sim_config.yaml` to `/tmp/sim_config.yaml` when the current
        working directory is `/tmp`).
    '''
    assert types.is_str_nonempty(pathname),\
        types.assert_nonstr_nonempty(pathname, 'Pathname')
    return path.abspath(path.expanduser(pathname))

# ....................{ MOVERS                             }....................
def move(pathname_source: str, pathname_target: str) -> None:
    '''
    Move the passed source to target path.

    Such path will be moved in a manner maximally preserving metadata (e.g.,
    owner, group, permissions, times, extended file system attributes).
    Likewise, if such source path is a symbolic link, such link (rather than its
    transitive target) will be moved and hence preserved.

    If either the source path does not exist *or* the target path already
    exists, an exception will be raised.
    '''
    assert types.is_str_nonempty(pathname_source),\
        types.assert_nonstr_nonempty(pathname_source, 'Source pathname')
    assert types.is_str_nonempty(pathname_target),\
        types.assert_nonstr_nonempty(pathname_target, 'Target pathname')

    # Log such move in a contextual manner.
    loggers.log_debug(
        'Moving %s "%s" to "%s".',
        get_type_label(pathname_source), pathname_source, pathname_target)

    # Raise an exception unless the source path exists.
    die_unless_path(pathname_source)

    # Raise an exception if the target path already exists. This is essential to
    # sane cross-platform semantics. shutil.move() is implemented in terms of
    # os.rename(), which overwrites such path if such path is a non-directory
    # file depending on platform. In most cases, that's bad.
    die_if_path(pathname_target)

    # Perform such move.
    shutil.move(pathname_source, pathname_target)

# --------------------( WASTELANDS                         )--------------------
# def die_unless_dirname_empty(pathname: str) -> None:
#     '''
#     Raise an exception unless the passed pathname is a pure basename.
#
#     See Also
#     ----------
#     `is_dirname_empty()`
#         For further details.
#     '''
#     if not is_dirname_empty(pathname):
#         raise BetseExceptionPath(
#             'Path "{}" contains directory separators.'.format(pathname))

# def get_dirname_or_none(pathname: str) -> str:
#     '''
#     Get the *dirname* (i.e., parent directory) of the passed path if such path
#     has a dirname or None otherwise.
#     '''
#     return get_dirname_or_empty(pathname) or None

    # Strip the prefixing "." from such filetype if any.
    # filetype = strs.remove_prefix_if_any(filetype, '.')
    #
    # # Test such filetype.
    # return get_filetype(pathname) == filetype

    # Such filetype may be either prefixed by `.` *or* not prefixed by `.`, but
    # should otherwise contain *no* `.` characters (e.g., `.gz` is accepted but
    # `.tar.gz` is *not*). In either case, this function operates as expected.

    # While there are more efficient implementations, the simplest should be
    # fine... for now.
    # return pathname == path.get_basename(pathname)

    # assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    # assert len(pathname), 'Pathname empty.'
    #
    # # Such dirname.
    # dirname = path.dirname(pathname)
    #
    # # Get such dirname. Since path.dirname() returns the empty string rather
    # # than None for paths without a dirname, convert the former to the latter.
    # return dirname if dirname else None
#FUXME: Obviously insufficient. This should also test whether such path is a

# ....................{ REMOVERS                           }....................
# def remove(filename: str) -> None:
#     '''
#     Remove the passed file.
#
#     If such file does *not* exist, an exception is raised.
#     '''
#     assert isinstance(filename, str), '"{}" not a string.'.format(filename)
#     os.remove(filename)

# def __init__():
#     '''
#
#     To support caller-specific exception handling, this function does *not*
#     validate such constants (e.g., by creating non-existent directories). See
#     `paths.init()` for such functionality.
#     '''
#     # Declare such constants to be globals, permitting their modification below.
#     global DATA_DIRNAME, DOT_DIRNAME
