#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility and convenience functions for application-specific :mod:`setuptools`
subcommands.

Design
----------
This method intentionally duplicates existing utility functions provided by the
:mod:`betse.util` subpackage. While duplication breaks DRY ("Don't Repeat
Yourself") and hence is usually harmful, there are valid reasons to do so here.
Namely, :mod:`betse.util` functionality:

* Assumes BETSE to be available. While this is certainly the case when this
  file resides in the BETSE codebase, this is *not* necessarily the case when
  this file is copied into and hence resides in the codebases of other projects
  (e.g., BETSEE). In these projects, BETSE is merely yet another dependency
  that is typically unavailable at installation time.
* Raises BETSE-specific exceptions rooted at the BETSE-specific
  :class:`betse.exception.BetseException` superclass. :mod:`setuptools`
  subcommands, on the other hand, are expected to only raise
  :mod:`distutils`-specific exceptions rooted at the :mod:`distutils`-specific
  :class:`DistutilsError` superclass.
* Assumes logging to be configured. :mod:`setuptools`, however, assumes
  logging to *not* be configured -- and provides no assistance in doing so.
* Could theoretically import third-party dependencies unavailable at
  :mod:`setuptools subcommand time (e.g., due to the ``install`` or ``develop``
  subcommands *not* having been run yet). While no :mod:`betse.util` submodules
  should do so, the horrid possibility remains.

Since duplicating these functions here is no significant maintenance burden
*and* since attempting to reuse these functions here would introduce spurious
side effects, we adopt the former approach.
'''

# ....................{ IMPORTS                           }....................
import os, platform, shutil, subprocess, sys, time
from distutils.errors import DistutilsFileError
from os import path

# ....................{ EXCEPTIONS ~ path                 }....................
def die_unless_basename(pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path is a **basename** (i.e., contains
    no platform-specific directory separator characters).
    '''

    # If this path is not a basename, fail.
    if not is_basename(pathname):
        # If no such message was passed, default this message.
        if not exception_message:
            exception_message = (
                'Pathname "{}" not a basename '
                '(i.e., either empty or '
                'contains a directory separator).'.format(pathname))
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsFileError(exception_message)


def die_unless_dir_or_not_found(
    pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path is either an existing directory
    *or* does not exist (i.e., if this path is an existing non-directory).
    '''
    # If such path is an existing non-directory, fail.
    if is_path(pathname) and not is_dir(pathname):
        # If no such message was passed, default such message.
        if not exception_message:
            if is_file(pathname):
                exception_message =\
                    'Directory "{}" already an existing file.'.format(pathname)
            elif is_symlink(pathname):
                exception_message =\
                    'Directory "{}" already an existing symbolic link.'.format(
                        pathname)
            else:
                exception_message = 'Path "{}" not a directory.'.format(
                    pathname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise this exception.
        raise DistutilsFileError(exception_message)


def die_unless_file_or_not_found(
    pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path is either an existing non-special
    file _or_ does not exist (e.g., if such path is an existing directory).
    '''

    # If this path exists and is *NOT* an existing non-special file, fail.
    if is_path(pathname) and not is_file(pathname):
        # If no such message was passed, default this message.
        if not exception_message:
            if is_dir(pathname):
                exception_message = (
                    'File "{}" already an existing directory.'.format(
                        pathname))
            elif is_symlink(pathname):
                exception_message = (
                    'File "{}" already an existing symbolic link.'.format(
                        pathname))
            else:
                exception_message = 'Path "{}" not a file.'.format(pathname)
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsFileError(exception_message)


def die_unless_path(pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path exists.
    '''

    # If this path is not found, fail.
    if not is_path(pathname):
        # If no such message was passed, default this message.
        if not exception_message:
            exception_message = 'Path "{}" not found.'.format(pathname)
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsFileError(exception_message)


def die_unless_dir(dirname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed directory exists.
    '''

    # If this directory is not found, fail.
    if not is_dir(dirname):
        # If no such message was passed, default such message.
        if not exception_message:
            exception_message = 'Directory "{}" not found.'.format(dirname)
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception. Since there exists no
        # "DistutilsDirError" class, the next best thing is raised.
        raise DistutilsFileError(exception_message)


def die_unless_file(filename: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed non-special file exists.
    '''

    # If such file is not found, fail.
    if not is_file(filename):
        # If no such message was passed, default such message.
        if not exception_message:
            exception_message = 'File "{}" not found.'.format(filename)
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsFileError(exception_message)


def die_unless_symlink(filename: str) -> None:
    '''
    Raise an exception unless the passed symbolic link exists.
    '''
    assert isinstance(filename, str), (
        '"{}" not a string.'.format(filename))

    if not is_symlink(filename):
        raise DistutilsFileError(
            'Symbolic link "{}" not found.'.format(filename))

# ....................{ TESTERS ~ os                      }....................
def is_os_posix() -> bool:
    '''
    `True` if the current operating system complies with POSIX standards (e.g.,
    as required for POSIX-compliant symbolic link support).

    Typically, this implies this system to _not_ be vanilla Microsoft Windows
    (i.e., to be either a Cygwin-enabled Windows terminal *or* a genuine
    POSIX-compliant system).
    '''
    return os.name == 'posix'
    # return False


def is_os_os_x() -> bool:
    '''
    `True` if the current operating system is Apple OS X.
    '''
    return platform.system() == 'Darwin'

# ....................{ TESTERS ~ os : windows            }....................
def is_os_windows() -> bool:
    '''
    `True` if the current operating system is Microsoft Windows.

    This function reports `True` for both vanilla and Cygwin Microsoft Windows.
    '''
    return is_os_windows_vanilla() or is_os_windows_cygwin()


def is_os_windows_cygwin() -> bool:
    '''
    `True` if the current operating system is **Cygwin Microsoft Windows**
    (i.e., running the Cygwin POSIX compatibility layer).
    '''
    return sys.platform == 'cygwin'


def is_os_windows_vanilla() -> bool:
    '''
    `True` if the current operating system is **vanilla Microsoft Windows**
    (i.e., _not_ running the Cygwin POSIX compatibility layer).
    '''
    return sys.platform == 'win32'

# ....................{ TESTERS ~ path                    }....................
def is_basename(pathname: str) -> bool:
    '''
    `True` only if the passed path is a **basename** (i.e., is a non-empty
    string containing no platform-specific directory separator characters).
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)

    return len(pathname) and pathname != path.basename(pathname)


def is_path(pathname: str) -> bool:
    '''
    `True` only if the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)

    return path.exists(pathname)


def is_dir(pathname: str) -> bool:
    '''
    `True` only if the passed directory exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)

    return path.isdir(pathname)


def is_file(pathname: str) -> bool:
    '''
    `True` only if the passed path is an existing non-directory file exists
    *after* following symbolic links.

    Versus `path.isfile()`
    ----------
    This function intrinsically differs from the standard `path.isfile()`
    function. While the latter returns `True` only for non-special files and
    hence `False` for all non-directory special files (e.g., device nodes,
    sockets), this function returns `True` for *all* non-directory files
    regardless of whether such files are special or not.

    **Why?** Because this function complies with POSIX semantics, whereas
    `path.isfile()` does *not*. The specialness of non-directory files is
    usually irrelevant; in general, it only matters whether such files are
    directories or not. For example, the external command `rm` removes only
    non-directory files (regardless of specialness) while the external command
    `rmdir` removes only empty directories.
    '''
    return is_path(pathname) and not is_dir(pathname)


def is_symlink(filename: str) -> bool:
    '''
    `True` if the passed symbolic link exists.

    `False` is returned if the passed symbolic link exists but the current user
    has insufficient privelages to follow such link. This may constitute a bug
    in the underlying `path.islink()` function.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    return path.islink(filename)

# ....................{ GETTERS ~ io : file               }....................
def get_chars(filename: str, encoding: str = 'utf-8') -> str:
    '''
    String of all characters contained in the plaintext file with the passed
    filename encoded with the passed encoding.

    Parameters
    ----------
    filename : str
        Relative or absolute path of the plaintext text to be read.
    encoding : optional[str]
        Name of the encoding to be used. Defaults to UTF-8.

    Returns
    ----------
    str
        String of all characters decoded from this file's byte content.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    assert isinstance(encoding, str), '"{}" not a string.'.format(encoding)

    with open(filename, mode='rt', encoding=encoding) as text_file:
        return text_file.read()

# ....................{ GETTERS ~ metadata                }....................
def get_description() -> str:
    '''
    Human-readable multiline description of this application in
    reStructuredText (reST) format.

    To minimize synchronization woes, this description is identical to the
    contents of the :doc:`/README.rst` file. When submitting this application
    package to PyPI, this description is re-used verbatim as this package's
    front matter.

    Caveats
    ----------
    This function is I/O intensive and hence should be called sparingly --
    ideally, only once by this application's top-level ``setup.py`` script.
    '''

    # Relative path of this application's front-facing documentation in
    # reStructuredText format, required by PyPI. This path resides outside this
    # application's package tree and hence is inlined here rather than provided
    # by the "betsee.guimetaapp" submodule.
    DESCRIPTION_FILENAME = 'README.rst'

    # Description read from this description file.
    try:
        description = get_chars(DESCRIPTION_FILENAME)
        # print('description: {}'.format(_DESCRIPTION))
    # If this file is *NOT* readable, print a non-fatal warning and reduce this
    # description to the empty string. While unfortunate, this description is
    # *NOT* required for most operations and hence mostly ignorable.
    except Exception as exception:
        description = ''
        output_warning(
            'Description file "{}" not found or not readable:\n{}'.format(
                DESCRIPTION_FILENAME, exception))

    # Retcurn this description.
    return description

# ....................{ GETTERS ~ path                    }....................
def get_path_canonicalized(pathname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed path.

    Specifically (in order):

    * Perform **tilde expansion,** replacing a `~` character prefixing such
      path by the absolute path of the current user's home directory.
    * Perform **path normalization,** thus:
      * Collapsing redundant separators (e.g., converting `//` to `/`).
      * Converting relative to absolute path components (e.g., converting `../`
        to the name of the parent directory of such component).
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.abspath(path.expanduser(pathname))


def get_path_dirname(pathname: str) -> str:
    '''
    Get the **dirname** (i.e., parent directory) of the passed path if such
    path has a dirname or raise an exception otherwise.
    '''
    # Get such dirname. Since the path.dirname() function returns the empty
    # string for paths containing no directory separators and hence having no
    # dirnames, assert such return value to be non-empty.
    dirname = path.dirname(pathname)
    assert len(dirname), 'Pathname "{}" dirname empty.'.format(pathname)
    return dirname

# ....................{ GETTERS ~ path : filetype         }....................
def get_path_filetype(pathname: str) -> str:
    '''
    Get the **last filetype** (i.e., last `.`-prefixed substring of the
    basename *not* including such `.`) of the passed path if this path has a
    filetype _or_ `None` otherwise.

    If this path contains multiple filetypes (e.g., `odium.reigns.tar.gz`),
    this function returns only the last filetype.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Such filetype. (Yes, splitext() is exceedingly poorly named.)
    filetype = path.splitext(pathname)[1]

    # Get such filetype, stripping the prefixing "." from the string returned
    # by the prior call if such path has a filetype or returning None
    # otherwise.
    return filetype[1:] if filetype else None

# ....................{ SANITIZERS ~ metadata             }....................
def sanitize_classifiers(
    classifiers: list,
    python_version_min_parts: tuple,
    python_version_minor_max: int,
) -> list:
    '''
    List of all PyPI-specific trove classifier strings synopsizing this
    application, manufactured by appending classifiers synopsizing this
    application's support for Python major versions (e.g.,
    ``Programming Language :: Python :: 3.6``, a classifier implying this
    application to successfully run under Python 3.6) to the passed list.

    Parameters
    ----------
    classifiers : list
        List of all PyPI-specific trove classifier strings to be sanitized.
    python_version_min_parts : tuple
        Minimum fully-specified version of Python required by this application
        as a tuple of integers (e.g., ``(3, 5, 0)`` if this application
        requires at least Python 3.5.0).
    python_version_minor_max : int
        Maximum minor stable version of the current Python 3.x mainline (e.g.,
        ``9`` if Python 3.9 is the most recent stable version of Python 3.x).

    Returns
    ----------
    list
        List of all sanitized PyPI-specific trove classifier strings.
    '''
    assert isinstance(classifiers, list), '"{}" not a list.'.format(
        classifiers)
    assert isinstance(python_version_min_parts, tuple), (
        '"{}" not a tuple.'.format(python_version_min_parts))
    assert isinstance(python_version_minor_max, int), (
        '"{}" not an integer.'.format(python_version_minor_max))

    # Major version of Python required by this application.
    PYTHON_VERSION_MAJOR = python_version_min_parts[0]

    # List of classifiers to return, copied from the passed list for safety.
    classifiers_sane = classifiers[:]

    # For each minor version of Python 3.x supported by this application,
    # formally classify this version as such.
    for python_version_minor in range(
        python_version_min_parts[1], python_version_minor_max):
        classifiers.append(
            'Programming Language :: Python :: {}.{}'.format(
                PYTHON_VERSION_MAJOR, python_version_minor,))
    # print('classifiers: {}'.format(_CLASSIFIERS))

    # Return this sanitized list of classifiers.
    return classifiers_sane

# ....................{ SANITIZERS ~ path                 }....................
def sanitize_command_basename(command_basename: str) -> str:
    '''
    Convert the passed platform-agnostic command basename (e.g., ``pytest``)
    into a platform-specific command basename (e.g., ``pytest.exe``).

    If the passed basename contains a directory separator and hence is *not* a
    basename, an exception is raised. Else, under:

    * Windows, the passed basename is appended by ``.exe``. To avoid confusion
      with non-Windows executables in the current ``${PATH}`` when running under
      Wine emulation, only Windows executables are accepted when running under
      Windows.
    * All other platforms, the passed basename is returned as is.
    '''
    assert isinstance(command_basename, str), (
        '"{}" not a string.'.format(command_basename))
    assert len(command_basename), 'Command basename empty.'

    # If this pathname is *NOT* a basename, raise an exception.
    die_unless_basename(command_basename)

    # If this is Windows *AND* this basename has no filetype, suffix this
    # basename by ".exe".
    if is_os_windows() and get_path_filetype(command_basename) is None:
        # print('command basename "{}" filetype; {}'.format(
        #     command_basename, get_path_filetype(command_basename)))
        return command_basename + '.exe'

    # Else, return this basename as is.
    return command_basename

# ....................{ OUTPUTTERS                        }....................
def output_warning(*warnings) -> None:
    '''
    Print the passed warning messages to standard error.
    '''
    print('WARNING: ', *warnings, file=sys.stderr)
