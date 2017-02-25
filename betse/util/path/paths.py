#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level path facilities specific to neither directories nor non-directory
files.

This module is named `paths` rather than `path` to avoid conflict with the stock
`path` module of the `os` package.
'''

# ....................{ IMPORTS                            }....................
import errno, os, shutil
from betse.exceptions import BetsePathException
from betse.util.io.log import logs
from betse.util.type.types import type_check, ContainerType
from os import path

# ....................{ CONSTANTS                          }....................
INVALID_PATHNAME = '\0'
'''
Pathname guaranteed to be invalid on all supported platforms.

This pathname consists of a single null byte. The operating system kernels for
both Windows _and_ all POSIX-compatible platforms (e.g., Linux, OS X) prohibit
pathnames containing null bytes by returning non-zero failure codes from
filesystem-centric kernel functions passed such pathnames. Python translates
these errors into raised exceptions. Under Linux, for example, a `TypeError`
exception of `"embedded NUL character"` is raised.
'''

# ....................{ EXCEPTIONS ~ path                  }....................
def die_if_path(*pathnames) -> None:
    '''
    Raise an exception if any passed path exists.
    '''

    for pathname in pathnames:
        if is_path(pathname):
            raise BetsePathException(
                'Path "{}" already exists.'.format(pathname))


def die_unless_path(*pathnames) -> None:
    '''
    Raise an exception unless all passed paths exist.
    '''

    for pathname in pathnames:
        if not is_path(pathname):
            raise BetsePathException(
                'Path "{}" not found or unreadable.'.format(pathname))

# ....................{ EXCEPTIONS ~ basename              }....................
def die_if_basename(pathname: str) -> None:
    '''
    Raise an exception if the passed path contains no directory separators.

    See Also
    ----------
    :func:`is_basename`
        Further details.
    '''

    if is_basename(pathname):
        raise BetsePathException(
            'Path "{}" contains no directory separators '
            "(i.e., '{}' characters).".format(pathname, os.sep))


def die_unless_basename(pathname: str) -> None:
    '''
    Raise an exception if the passed path contains one or more directory
    separators.

    See Also
    ----------
    is_basename()
        For further details.
    '''

    if not is_basename(pathname):
        raise BetsePathException(
            'Path "{}" contains one or more directory separators'
            "(i.e., '{}' characters).".format(pathname, os.sep))

# ....................{ TESTERS                            }....................
@type_check
def is_path(pathname: str) -> bool:
    '''
    `True` only if the passed path exists _without_ following this path if this
    path is an existing symbolic link.

    If this path is an existing **dangling symbolic link** (i.e., symbolic link
    whose target no longer exists, also commonly referred to as a broken
    symbolic link), this function still returns `True`.
    '''

    # Call path.lexists() rather than path.exists(), as the latter returns False
    # for dangling symbolic links -- which is entirely irrelevant for most
    # contexts. Under POSIX semantics, dangling symbolic links are essential to
    # common usage patterns and should *NOT* be discriminated against here.
    return path.lexists(pathname)

# ....................{ TESTERS ~ absolute                 }....................
@type_check
def is_absolute(pathname: str) -> bool:
    '''
    `True` only if the passed path is absolute.

    The definition of "absolute" depends on the current operating system. Under:

    * POSIX-compatible systems, absolute paths are prefixed by `/`.
    * Microsoft Windows, absolute paths are prefixed by an optional drive
      indicator (e.g., `C:`) followed by `\`.
    '''

    return path.isabs(pathname)


def is_relative(pathname: str) -> bool:
    '''
    `True` only if the passed path is relative.

    The definition of "relative" depends on the current operating system. Under:

    * POSIX-compatible systems, relative paths are _not_ prefixed by `/`.
    * Microsoft Windows, relative paths are _not_ prefixed by an optional drive
      indicator (e.g., `C:`) followed by `\`.
    '''
    return not is_absolute(pathname)

# ....................{ TESTERS ~ pathname                 }....................
@type_check
def is_basename(pathname: str) -> bool:
    '''
    `True` only if the passed pathname is a **pure basename** (i.e., contains no
    directory separators and hence no directory components).
    '''

    return path.sep not in pathname


@type_check
def is_pathname(pathname: str) -> bool:
    '''
    `True` only if the passed string is a valid pathname (either absolute or
    relative) for the root filesystem of the current platform.

    Under:

    * POSIX-compatible OSes:
      * Valid pathnames are non-empty strings containing:
        * No null byte characters.
        * No `/`-delimited path component longer than 255 characters.
      * The root filesystem is the filesystem mounted to the root directory `/`.
    * Microsoft OSes:
      * Valid pathnames are non-empty strings satisfying `various constraints
        <https://msdn.microsoft.com/en-us/library/windows/desktop/aa365247%28v=vs.85%29.aspx>`_
        too numerous to document here.
      * The root filesystem is the filesystem to which this instance of Windows
        was installed, also given by the `%HOMEDRIVE%` environment variable.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.brand import windows
    from betse.util.os.brand.windows import WindowsErrorType

    # If this is the empty string, it cannot by definition be a valid pathname.
    if not pathname:
        return False

    # The only cross-platform and -filesystem portable means of validating
    # pathnames is to parse exceptions raised by the kernel-dependent os.stat()
    # or os.lstat() functions for metadata indicating invalid pathname syntax.
    # All other approaches (e.g., regex string matching) fail for common edge
    # cases. See also:
    #     https://stackoverflow.com/a/34102855/2809027
    try:
        # Strip this pathname's Windows drive specifier (e.g., "C:\") if any.
        # Since Windows prohibits path components from containing ":"
        # characters, failing to strip this ":"-suffixed prefix would
        # erroneously invalidate all valid absolute Windows pathnames.
        _, pathname = path.splitdrive(pathname)

        # Absolute path of a directory guaranteed to exist.
        #
        # To avoid race conditions with external processes concurrently
        # modifying the filesystem, the passed pathname cannot be tested as is.
        # Only path components split from this pathname are safely testable.
        # Why? Because os.stat() and os.lstat() raise "FileNotFoundError"
        # exceptions when passed pathnames residing in non-existing directories
        # regardless of whether these pathnames are invalid or not. Directory
        # existence takes precedence over pathname invalidity. Hence, the only
        # means of testing whether pathnames are invalid or not is to:
        #
        # 1. Split the passed pathname into path components (e.g.,
        #    "/foo/bar" into "['', 'foo', 'bar']").
        # 2. For each path component:
        #    1. Join the pathname of a directory guaranteed to exist and the
        #       current path component into a new pathname (e.g., "/bar").
        #    2. Pass that pathname to os.stat() or os.lstat(). If that
        #       pathname and hence current path component is invalid, this
        #       call is guaranteed to raise an exception exposing the type
        #       of invalidity rather than a generic "FileNotFoundError"
        #       exception. Why? Because that pathname resides in an
        #       existing directory. Circular logic is circular.
        #
        # Is there a directory guaranteed to exist? There is, but typically only
        # one: the root directory for the current filesystem. Passing pathnames
        # residing in any other directories to os.stat() or os.lstat() invites
        # mind-flaying race conditions, even for directories previously tested
        # to exist. External processes cannot be prevented from concurrently
        # removing those directories after those tests have been performed but
        # before those pathnames are passed to os.stat() or os.lstat().
        #
        # Did we mention this should be shipped with Python already?
        root_dirname = get_root_dirname()
        assert os.path.isdir(root_dirname)   # ...Murphy and her dumb Law

        # Test whether each path component split from this pathname is valid
        # or not. Most path components will *NOT* actually physically exist.
        for pathname_part in pathname.split(path.sep):
            try:
                os.lstat(root_dirname + pathname_part)
            # If an OS-specific exception is raised, its error code indicates
            # whether this pathname is valid or not. Unless this is the case,
            # this exception implies an ignorable kernel or filesystem complaint
            # (e.g., path not found or inaccessible).
            #
            # Only the following exceptions indicate invalid pathnames:
            #
            # * Instances of the Windows-specific "WindowsError" class
            #   defining the "winerror" attribute whose value is
            #   "ERROR_INVALID_NAME". Under Windows, "winerror" is more
            #   fine-grained and hence useful than the generic "errno"
            #   attribute. When a too-long pathname is passed, for example,
            #   "errno" is "ENOENT" (i.e., no such file or directory) rather
            #   than "ENAMETOOLONG" (i.e., file name too long).
            # * Instances of the cross-platform "OSError" class defining the
            #   generic "errno" attribute whose value is either:
            #   * Under most POSIX-compatible OSes, "ENAMETOOLONG".
            #   * Under some edge-case OSes (e.g., SunOS, *BSD), "ERANGE".
            except OSError as exc:
                # True only if this pathname is invalid (as detailed above).
                is_pathname_invalid = (
                    windows.is_exception_pathname_invalid(exc)
                    if isinstance(exc, WindowsErrorType)
                    else exc.errno in {errno.ENAMETOOLONG, errno.ERANGE})

                # If this pathname is invalid, log a warning and return False.
                if is_pathname_invalid:
                    logs.log_warning(
                        'Pathname "{}" invalid: {}'.format(
                            pathname, exc.strerror))
                    return False
    # If a "TypeError" exception was raised, it almost certainly has the
    # error message "embedded NUL character" indicating an invalid pathname.
    except TypeError as exc:
        logs.log_warning(
            'Pathname "{}" invalid: {}'.format(pathname, str(exc)))
        return False
    # If no exception was raised, all path components and hence this
    # pathname itself are valid. (Praise be to the curmudgeonly python.)
    else:
        return True
    # If any other exception was raised, this is an unrelated fatal issue
    # (e.g., a bug). Permit this exception to unwind the call stack.

# ....................{ TESTERS ~ filetype                 }....................
@type_check
def is_filetype_equals(pathname: str, filetype: str) -> bool:
    '''
    `True` only if the passed pathname has a filetype _and_ the **last
    filetype** (i.e., last `.`-prefixed substring of the basename) of this
    pathname is the passed filetype.

    The passed filetype may contain arbitrarily many `.` characters, including
    an optional prefixing `.`. In all cases, this function behaves as expected.

    Parameters
    ----------
    pathname : str
        Pathname to be tested.
    filetype : str
        Filetype to test whether this pathname is suffixed by.

    Returns
    ----------
    bool
        `True` only if this pathname has this filetype.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import strs

    # Test this filetype, prefixed by "." unless already prefixed.
    return strs.is_suffix(
        pathname, strs.add_prefix_unless_found(filetype, '.'))


@type_check
def is_filetype_undotted_in(
    pathname: str, filetypes_undotted: ContainerType) -> bool:
    '''
    `True` only if the passed pathname has a filetype _and_ the **last undotted
    filetype** (i.e., last `.`-prefixed substring of the basename excluding this
    `.`) of this pathname is an element of the passed container.

    Parameters
    ----------
    pathname : str
        Pathname to be tested.
    filetypes: ContainerType
        Container of all undotted filetypes to test this pathname against. Each
        element of this container should be a string _not_ prefixed the `.`
        character.

    Returns
    ----------
    bool
        `True` only if this pathname has a filetype _and_ the last filetype of
        this pathname is an element of this container.
    '''

    # Undotted filetype of this filename if any or None otherwise
    filetype_undotted = get_filetype_undotted_or_none(pathname)

    # Return True only if...
    return (
        # This pathname is suffixed by a filetype.
        filetype_undotted is not None and
        # This filetype is in this container.
        filetype_undotted in filetypes_undotted
    )

# ....................{ GETTERS                            }....................
@type_check
def get_basename(pathname: str) -> str:
    '''
    **Basename** (i.e., last component) of the passed path.
    '''

    return path.basename(pathname)


def get_type_label(pathname: str) -> str:
    '''
    Get a lowercase human-readable label describing the type of the passed
    existing path.

    For brevity, this label comprises one lowercase noun optionally preceded by
    a lowercase adjective unambiguously identifying such path's type:

    * `symbolic link` for an existing symbolic link, regardless of whether that
      link is **dangling** (i.e., refers to a non-existent path) or not.
    * `directory` for an existing directory.
    * `special file` for all other **special files** (e.g., device node, named
      pipe, socket).
    * `file` for an existing regular (i.e., non-special) file).

    If this path does _not_ exist, an exception is raised.
    '''

    # Avoid circular import dependencies.
    from betse.util.paths import dirs, files

    # Raise an exception unless this path exists.
    die_unless_path(pathname)

    # Get this type.
    if files.is_file(pathname):
        if files.is_symlink(pathname):
            return 'symbolic link'
        elif files.is_special(pathname):
            return 'special file'
        else:
            return 'file'
    elif dirs.is_dir(pathname):
        return 'directory'
    # Else, this existing path is neither a file or directory and is hence a
    # shambling thing that should not be. Panic, for space-time continuity
    # itself is now in question.
    else:
        return BetsePathException(
            'Path "{}" type unrecognized.'.format(pathname))

    # This should never happen. That means it will in several seconds.
    return None

# ....................{ GETTERS ~ dirname                  }....................
def get_root_dirname() -> str:
    '''
    Absolute path of the root directory suffixed by a directory separator.

    The definition of "root directory" conditionally depends on the current
    platform. If this platform is:

    * POSIX-compatible (e.g., Linux, OS X), this is simply `/`.
    * Microsoft Windows, this is the value of the `%HOMEDRIVE%` environment
      variable. This is the `:`-suffixed letter of the drive to which Windows
      was originally installed -- typically but _not_ necessarily `C:\`.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # Return this path.
    if oses.is_windows_vanilla():
        return os.environ.get('HOMEDRIVE', 'C:') + os.path.sep
    else:
        return path.sep


def get_dirname(pathname: str) -> str:
    '''
    **Dirname** (i.e., parent directory) of the passed path if this path has a
    dirname *or* raise an exception otherwise.
    '''

    # Raise an exception unless this path has a dirname.
    die_if_basename(pathname)

    # This dirname.
    dirname = get_dirname_or_empty(pathname)

    # Assert this dirname's non-emptiness. Technically, the above call *SHOULD*
    # have ensured this. This is a sufficiently critical function, however, to
    # warrant asserting this fact.
    assert len(dirname), 'Pathname "{}" dirname empty.'.format(pathname)
    return dirname


def get_dirname_or_current_dirname(pathname: str) -> str:
    '''
    **Dirname** (i.e., parent directory) of the passed path if this path has a
    dirname *or* the current working directory otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs
    dirname = get_dirname_or_empty(pathname)
    return dirname if dirname else dirs.get_current_dirname()


@type_check
def get_dirname_or_empty(pathname: str) -> str:
    '''
    **Dirname** (i.e., parent directory) of the passed path if this path has a
    dirname *or* the empty string otherwise.
    '''

    return path.dirname(pathname)

# ....................{ GETTERS ~ filetype                 }....................
@type_check
def get_pathname_sans_filetypes(pathname: str) -> str:
    '''
    Passed path without all suffixing ``.``-prefixed filetypes (including these
    prefixes) if this path has any filetypes *or* this path as is otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import strs

    # Dirname of this path if any or the empty string otherwise.
    dirname = get_dirname_or_empty(pathname)

    # Basename of this path.
    basename = get_basename(pathname)

    # Strip all characters following the first "." from this basename.
    basename = strs.remove_suffix_with_prefix(basename, '.')

    # If this path contains a dirname, return the concatenation of this dirname
    # as is by this stripped basename.
    if dirname:
        return path.join(dirname, basename)
    # Else, return merely this stripped basename.
    else:
        return basename


@type_check
def get_pathname_sans_filetype(pathname: str) -> str:
    '''
    Passed path without the last ``.``-prefixed filetype (including this prefix)
    if this path has a filetype *or* this path as is otherwise.
    '''

    return path.splitext(pathname)[0]

# ....................{ GETTERS ~ filetype : dotted        }....................
@type_check
def get_filetype_dotted_or_none(pathname: str) -> str:
    '''
    ``.``-prefixed **last filetype** (i.e., last ``.``-prefixed substring of the
    basename) of the passed pathname if this pathname is suffixed by a filetype
    *or* ``None`` otherwise.

    If this pathname is suffixed by two or more filetypes (e.g.,
    ``odium.reigns.tar.gz``), only the last such filetype is returned.
    '''

    # "."-prefixed filetype of this pathname.
    filetype = path.splitext(pathname)[1]

    # Return this string as is if non-empty or "None" otherwise.
    return filetype or None

# ....................{ GETTERS ~ filetype : undotted      }....................
@type_check
def get_filetype_undotted(pathname: str) -> str:
    '''
    **Last filetype** (i.e., last ``.``-prefixed substring of the basename
    excluding this ``.``) of the passed pathname.

    Raises
    ----------
    :exc:`BetsePathException`
        If this pathname is *not* suffixed by a filetype.

    See Also
    ----------
    :func:`get_filetype_dotted_or_none`
        Further details.
    '''

    # Filetype of this pathname if any or "None" otherwise.
    filetype = get_filetype_undotted_or_none(pathname)

    # If this pathname has no filetype, raise an exception.
    if filetype is None:
        raise BetsePathException('Path "{}" has no filetype.'.format(pathname))

    # Else, return this filetype.
    return filetype


@type_check
def get_filetype_undotted_or_none(pathname: str) -> str:
    '''
    **Last filetype** (i.e., last ``.``-prefixed substring of the basename
    excluding this ``.``) of the passed pathname if this pathname is suffixed by
    a filetype *or* ``None`` otherwise.

    See Also
    ----------
    :func:`get_filetype_dotted_or_none`
        Further details.
    '''

    # "."-prefixed filetype of this pathname.
    filetype = get_filetype_dotted_or_none(pathname)

    # Return this string stripped of this prefix if non-empty or "None"
    # otherwise.
    return filetype[1:] if filetype is not None else None

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
    function _without_ adding functionality to that function -- principally to
    unify and hence simplify `import` statements in other modules.
    '''

    return path.join(*pathnames)

# ....................{ CANONICALIZERS                     }....................
@type_check
def canonicalize(pathname: str) -> str:
    '''
    **Canonical form** (i.e., unique absolute pathname *after* transitively
    resolving all symbolic links) of the passed path.

    Specifically (in order):

    . Transitively resolve all symbolic links, producing a pathname that either
      does not exist _or_ does exist but is not a symbolic link.
    . Perform **tilde expansion,** replacing a ``~`` character prefixing this
      path by the absolute path of the current user's home directory.
    . Perform **path normalization,** thus (in no particular order):
      * Collapsing redundant separators (e.g., converting ``//`` to ``/``).
      * Converting explicit relative to absolute path components (e.g.,
        converting ``../`` to the name of the parent directory of that
        component).
      * Converting implicit relative basenames to absolute pathnames (e.g.,
        converting ``sim_config.yaml`` to ``/tmp/sim_config.yaml`` when the
        current working directory is ``/tmp``).
    '''

    return path.realpath(path.expanduser(pathname))

# ....................{ MOVERS                             }....................
@type_check
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

    # Log such move in a contextual manner.
    logs.log_debug(
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
