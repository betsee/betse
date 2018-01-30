#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level pathname (e.g., basename, dirname, filetype) facilities.
'''

# ....................{ IMPORTS                            }....................
import errno, os
from betse.exceptions import BetsePathnameException
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, ContainerType, ModuleType, StrOrNoneTypes)
from os import path as os_path

# ....................{ CONSTANTS                          }....................
INVALID_PATHNAME = '\0'
'''
Pathname guaranteed to be invalid on all supported platforms.

This pathname consists of a single null byte. The operating system kernels for
both Windows _and_ all POSIX-compatible platforms (e.g., Linux, OS X) prohibit
pathnames containing null bytes by returning non-zero failure codes from
filesystem-centric kernel functions passed such pathnames. Python translates
these errors into raised exceptions. Under Linux, for example, a
:class:`TypeError` exception of message ``"embedded NUL character"`` is raised.
'''

# ....................{ EXCEPTIONS ~ absolute              }....................
def die_if_absolute(*pathnames: str) -> None:
    '''
    Raise an exception if any passed path is absolute (i.e., unless all passed
    paths are relative).

    See Also
    ----------
    :func:`is_absolute`
        Further details.
    '''

    for pathname in pathnames:
        if is_absolute(pathname):
            raise BetsePathnameException(
                'Pathname "{}" absolute rather than relative.'.format(pathname))


def die_if_relative(*pathnames: str) -> None:
    '''
    Raise an exception if any passed path is relative (i.e., unless all passed
    paths are absolute).

    See Also
    ----------
    :func:`is_relative`
        Further details.
    '''

    for pathname in pathnames:
        if is_relative(pathname):
            raise BetsePathnameException(
                'Pathname "{}" relative rather than absolute.'.format(pathname))

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
        raise BetsePathnameException(
            'Pathname "{}" contains no directory separators '
            "(i.e., '{}' characters).".format(pathname, os_path.sep))


def die_unless_basename(pathname: str) -> None:
    '''
    Raise an exception if the passed path contains one or more directory
    separators.

    See Also
    ----------
    :func:`is_basename`)
        Further details.
    '''

    if not is_basename(pathname):
        raise BetsePathnameException(
            'Pathname "{}" contains one or more directory separators'
            "(i.e., '{}' characters).".format(pathname, os_path.sep))

# ....................{ EXCEPTIONS ~ parent                }....................
def die_unless_parent(parent_dirname: str, child_pathname: str) -> bool:
    '''
    Raise an exception unless the second passed pathname is a child of (i.e.,
    prefixed by) the first passed dirname.

    Parameters
    ----------
    parent_dirname : str
        Absolute or relative pathname of the candidate parent directory.
    child_pathname : str
        Absolute or relative pathname of the candidate child path.

    See Also
    ----------
    :func:`is_parent`
        Further details.
    '''

    if not is_parent(parent_dirname, child_pathname):
        raise BetsePathnameException(
            'Pathname "{}" not in dirname "{}".'.format(
                child_pathname, parent_dirname,))

# ....................{ EXCEPTIONS ~ filetype              }....................
def die_unless_filetype_equals(pathname: str, filetype: str) -> None:
    '''
    Raise an exception unless the passed path has the passed filetype.

    See Also
    ----------
    :func:`is_filetype_equals`
        Further details.
    '''

    if not is_filetype_equals(pathname, filetype):
        raise BetsePathnameException(
            'Pathname "{}" filetype not "{}".'.format(pathname, filetype))

# ....................{ TESTERS                            }....................
@type_check
def is_pathname(pathname: str) -> bool:
    '''
    ``True`` only if the passed string is a valid pathname (either absolute or
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
    from betse import pathtree
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
        _, pathname = os_path.splitdrive(pathname)

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
        root_dirname = pathtree.get_root_dirname()
        assert os_path.isdir(root_dirname)   # ...Murphy and her dumb Law

        # Test whether each path component split from this pathname is valid
        # or not. Most path components will *NOT* actually physically exist.
        for pathname_part in pathname.split(os_path.sep):
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

# ....................{ TESTERS ~ absolute                 }....................
@type_check
def is_absolute(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is absolute.

    The definition of "absolute" depends on the current operating system. Under:

    * POSIX-compatible systems, absolute paths are prefixed by ``/``.
    * Microsoft Windows, absolute paths are prefixed by an optional drive
      indicator (e.g., ``C:``) followed by ``\``.
    '''

    return os_path.isabs(pathname)


def is_relative(pathname: str) -> bool:
    '''
    ``True`` only if the passed path is relative.

    The definition of "relative" depends on the current operating system. Under:

    * POSIX-compatible systems, relative paths are *not* prefixed by ``/``.
    * Microsoft Windows, relative paths are *not* prefixed by an optional drive
      indicator (e.g., ``C:``) followed by ``\``.
    '''

    return not is_absolute(pathname)

# ....................{ TESTERS ~ basename                 }....................
@type_check
def is_basename(pathname: str) -> bool:
    '''
    ``True`` only if the passed pathname is a **pure basename** (i.e., contains
    no directory separators and hence no directory components).
    '''

    return os_path.sep not in pathname

# ....................{ TESTERS ~ parent                   }....................
@type_check
def is_parent(parent_dirname: str, child_pathname: str) -> bool:
    '''
    ``True`` only if the second passed pathname is a child of (i.e., prefixed
    by) the first passed dirname.

    Parameters
    ----------
    parent_dirname : str
        Absolute or relative pathname of the candidate parent directory.
    child_pathname : str
        Absolute or relative pathname of the candidate child path.

    Returns
    ----------
    bool
        ``True`` only if this child path is a child of this parent directory.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Absolute or relative pathname of the candidate parent directory, suffixed
    # by this platform's directory separator to ensure this child is actually a
    # contained child of this parent directory.
    parent_dirname_suffixed = strs.add_suffix_unless_found(
        text=parent_dirname, suffix=os_path.sep)

    # Return true only if this child pathname is suffixed by this dirname.
    return strs.is_prefix(text=child_pathname, prefix=parent_dirname_suffixed)

# ....................{ TESTERS ~ filetype                 }....................
@type_check
def is_filetype_equals(pathname: str, filetype: str) -> bool:
    '''
    ``True`` only if the passed pathname has a filetype *and* the **last
    filetype** (i.e., last ``.``-prefixed substring of the basename) of this
    pathname is the passed filetype.

    The passed filetype may contain arbitrarily many ``.`` characters, including
    an optional prefixing ``.``. In all cases, this function behaves as
    expected.

    Parameters
    ----------
    pathname : str
        Pathname to be tested.
    filetype : str
        Filetype to test whether this pathname is suffixed by.

    Returns
    ----------
    bool
        ``True`` only if this pathname has this filetype.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Prefix this filetype by "." if needed.
    filetype = dot_filetype(filetype)

    # Return true only if this pathname is suffixed by this filetype.
    return strs.is_suffix(pathname, filetype)


@type_check
def is_filetype_undotted_in(
    pathname: str, filetypes_undotted: ContainerType) -> bool:
    '''
    ``True`` only if the passed pathname has a filetype *and* the **last
    undotted filetype** (i.e., last ``.``-prefixed substring of the basename
    excluding this ``.``) of this pathname is an element of the passed
    container.

    Parameters
    ----------
    pathname : str
        Pathname to be tested.
    filetypes: ContainerType
        Container of all undotted filetypes to test this pathname against. Each
        element of this container should be a string *not* prefixed the ``.``
        character.

    Returns
    ----------
    bool
        ``True`` only if this pathname has a filetype *and* the last filetype of
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

    return os_path.basename(pathname)

# ....................{ GETTERS ~ app                      }....................
@type_check
def get_app_pathname(package: ModuleType, pathname: str) -> str:
    '''
    Absolute path of the passed pathname relative to the absolute path of the
    directory providing the passed top-level package for this application.

    Specifically, this function returns:

    * If this application is a PyInstaller-frozen executable binary, the
      concatenation of (in order):

      #. The absolute path of the temporary directory containing all application
         data resources extracted from this binary by this executable's
         bootloader, as specified by the PyInstaller-specific private attribute
         ``_MEIPASS`` injected into the canonical :mod:`sys` module by the
         PyInstaller bootloader embedded in this binary. "And it's turtles all
         the way down."
      #. The passed relative pathname.

    * If this application is a :mod:`setuptools`-installed script wrapper, the
      result of querying :mod:`setuptools` for the absolute path of the passed
      relative pathname. In this case, this path will have been preserved as is
      in the :mod:`setuptools`-installed copy of this application in the package
      tree for the active Python interpreter.
    * Else, the concatenation of (in order):

      #. The absolute path of the directory providing this application package.
      #. The passed relative pathname.

      In this case, this application is typically either a
      :mod:`setuptools`-symlinked script wrapper *or* was invoked via the
      secretive ``python3 -m {package.__name__}`` command.

    Parameters
    ----------
    pathname : str
        Relative pathname of the path to be canonicalized.
    package : ModuleType
        Top-level package for this application.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import resources
    from betse.util.py import pyfreeze

    # If this application is frozen by PyInstaller, canonicalize this path
    # relative to the directory to which this application is unfrozen.
    if pyfreeze.is_frozen_pyinstaller():
        app_pathname = join(
            pyfreeze.get_app_dirname_pyinstaller(), pathname)
    # Else if this application is a setuptools-installed script wrapper,
    # canonicalize this path by deferring to the setuptools resource API.
    elif resources.is_dir(
        module_name=package.__name__, dirname=pathname):
        app_pathname = resources.get_pathname(
            module_name=package.__name__, pathname=pathname)
    # Else, the current application is either a setuptools-symlinked script
    # wrapper *OR* was invoked via the secretive "python3 -m betse" command. In
    # either case, this directory's path is directly obtainable relative to the
    # absolute path of the passed package.
    else:
        app_pathname = join(package.__path__, pathname)

    # Return this pathname.
    return app_pathname

# ....................{ GETTERS ~ dirname                  }....................
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


def get_dirname_or_cwd(pathname: str) -> str:
    '''
    **Dirname** (i.e., parent directory) of the passed path if this path has a
    dirname *or* the current working directory otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.shell import shelldir
    dirname = get_dirname_or_empty(pathname)
    return dirname if dirname else shelldir.get_cwd_dirname()


@type_check
def get_dirname_or_empty(pathname: str) -> str:
    '''
    **Dirname** (i.e., parent directory) of the passed path if this path has a
    dirname *or* the empty string otherwise.
    '''

    return os_path.dirname(pathname)

# ....................{ GETTERS ~ filetype                 }....................
def get_pathname_sans_filetypes(pathname: str) -> str:
    '''
    Passed path without all suffixing ``.``-prefixed filetypes (including these
    prefixes) if this path has any filetypes *or* this path as is otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Dirname of this path if any or the empty string otherwise.
    dirname = get_dirname_or_empty(pathname)

    # Basename of this path.
    basename = get_basename(pathname)

    # Strip all characters following the first "." from this basename.
    basename = strs.remove_suffix_with_prefix(basename, '.')

    # If this path contains a dirname, return the concatenation of this dirname
    # as is by this stripped basename.
    if dirname:
        return join(dirname, basename)
    # Else, return merely this stripped basename.
    else:
        return basename


@type_check
def get_pathname_sans_filetype(pathname: str) -> str:
    '''
    Passed path without the last ``.``-prefixed filetype (including this prefix)
    if this path has a filetype *or* this path as is otherwise.
    '''

    return os_path.splitext(pathname)[0]

# ....................{ GETTERS ~ filetype : undotted      }....................
@type_check
def get_filetype_dotted_or_none(pathname: str) -> StrOrNoneTypes:
    '''
    ``.``-prefixed **last filetype** (i.e., last ``.``-prefixed substring of the
    basename) of the passed pathname if this pathname is suffixed by a filetype
    *or* ``None`` otherwise.

    If this pathname is suffixed by two or more filetypes (e.g.,
    ``odium.reigns.tar.gz``), only the last such filetype is returned.
    '''

    # "."-prefixed filetype of this pathname.
    filetype = os_path.splitext(pathname)[1]

    # Return this string as is if non-empty or "None" otherwise.
    return filetype or None


def get_filetype_undotted(pathname: str) -> str:
    '''
    **Last filetype** (i.e., last ``.``-prefixed substring of the basename
    excluding this ``.``) of the passed pathname.

    Raises
    ----------
    :exc:`BetsePathnameException`
        If this pathname is suffixed by no filetype (i.e., if this pathname
        contains no ``.`` character).

    See Also
    ----------
    :func:`get_filetype_dotted_or_none`
        Further details.
    '''

    # Filetype of this pathname if any or "None" otherwise.
    filetype = get_filetype_undotted_or_none(pathname)

    # If this pathname has no filetype, raise an exception.
    if filetype is None:
        raise BetsePathnameException(
            'Path "{}" has no filetype.'.format(pathname))

    # Else, return this filetype.
    return filetype


@type_check
def get_filetype_undotted_or_none(pathname: str) -> StrOrNoneTypes:
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

# ....................{ DOTTERS                            }....................
@type_check
def dot_filetype(filetype: str) -> str:
    '''
    The passed filetype prefixed by a ``.`` character if needed (i.e., if not
    already prefixed by this character *or* this filetype as is otherwise).

    Parameters
    ----------
    filetype : str
        Filetype to be prefixed.

    Returns
    ----------
    str
        This filetype prefixed by a ``.`` character if *not* already prefixed.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Return this filetype prefixed by "." if needed.
    return strs.add_prefix_unless_found(text=filetype, prefix='.')


@type_check
def undot_filetype(filetype: str) -> str:
    '''
    The passed filetype prefixed by no ``.`` character.

    Parameters
    ----------
    filetype : str
        Filetype to be unprefixed.

    Returns
    ----------
    str
        Either:
        * If this filetype is prefixed by a ``.`` character, this filetype
          stripped of this prefix.
        * Else, this filetype unmodified.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Return this filetype prefixed by "." if needed.
    return strs.remove_prefix_if_found(text=filetype, prefix='.')

# ....................{ CANONICALIZERS                     }....................
@type_check
def join_and_canonicalize(*partnames: str) -> str:
    '''
    Canonical form of the path produced by joining (i.e., concatenating) the
    passed pathnames with the directory separator specific to the current
    platform.

    This higher-level function is a convenience wrapper encapsulating both the
    lower-level :func:`join` and :func:`canonicalize` functions.
    '''

    return canonicalize(join(*partnames))


@type_check
def canonicalize(pathname: str) -> str:
    '''
    **Canonical form** (i.e., unique absolute pathname *after* transitively
    resolving all symbolic links) of the passed path.

    Specifically, this function (in order):

    #. Transitively resolves all symbolic links, producing a pathname that
       #either:

       * Does not exist.
       * Does exist but is *not* a symbolic link.

    #. Performs **tilde expansion,** replacing a ``~`` character prefixing this
       path by the absolute path of the current user's home directory.

    #. Performs **path normalization,** which (in no particular order):

      * Collapses redundant separators (e.g., converting ``//`` to ``/``).
      * Converts explicit relative to absolute path components (e.g., converts
        ``../`` to the name of the parent directory of that component).
      * Converts implicit relative basenames to absolute pathnames (e.g.,
        converts ``sim_config.yaml`` to ``/tmp/sim_config.yaml`` when the
        current working directory is ``/tmp``).
    '''

    return os_path.realpath(os_path.expanduser(pathname))

# ....................{ JOINERS                            }....................
@type_check
def join(*partnames: str) -> str:
    '''
    Join (i.e., concatenate) the passed pathnames with the directory separator
    specific to the current platform.

    Caveats
    ----------
    This high-level wrapper should *always* be called in lieu of the low-level
    :func:`os.path.join` function, which sadly violates Python's core "Explicit
    is better than implicit." design principle. Specifically, to quote official
    :func:`os.path.join` documentation:

        If a component is an absolute path, all previous components are thrown
        away and joining continues from the absolute path component.

    Non-intuitive path mangling of that sort only promotes non-debuggable bugs.
    This function prevents these bugs by explicitly raising an exception if any
    passed pathname except the first is absolute.

    Parameters
    ----------
    partnames : tuple[str]
        Tuple of all pathnames to be joined. The first such pathname may be
        absolute or relative. All remaining pathnames *must* be relative.

    Returns
    ----------
    str
        Absolute or relative pathname joined from these pathnames.

    Raises
    ----------
    BetsePathnameException
        If any such pathname except the first is absolute.
    '''

    # If more than one partname was passed *AND* any such partname is absolute,
    # raise an exception. See above for commentary.
    if len(partnames) > 1:
        die_if_absolute(*partnames[1:])

    # Return the concatenation of these partnames.
    return os_path.join(*partnames)

# ....................{ RELATIVIZERS                       }....................
@type_check
def relativize(src_dirname: str, trg_pathname: str) -> str:
    '''
    Relative pathname of the second passed pathname relative to the first passed
    dirname.

    Parameters
    ----------
    src_dirname : str
        Absolute or relative pathname of the source directory.
    trg_pathname : str
        Absolute or relative pathname of the target path.

    Returns
    ----------
    str
        Absolute or relative pathname of the target path converted into a
        relative path relative to the source directory.
    '''

    return os_path.relpath(trg_pathname, start=src_dirname)
