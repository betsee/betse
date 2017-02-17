#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level directory facilities.
'''

# ....................{ IMPORTS                            }....................
import os, shutil
from betse.exceptions import BetseDirException
from betse.util.io.log import logs
from betse.util.type.types import type_check, GeneratorType
from contextlib import contextmanager
from os import path

# ....................{ GLOBALS                            }....................
# There exist only two possible directory separators for all modern platforms.
# Hence, this reliably suffices with no error handling required.
SEPARATOR_REGEX = r'/' if path.sep == '/' else r'\\'
'''
Regular expression matching the directory separator specific to the current
platform.

Specifically, under:

* Microsoft Windows, this is `\\\\`.
* All other platforms, this is `/`.
'''

# ....................{ EXCEPTIONS ~ unless                }....................
def die_unless_dir(*dirnames) -> None:
    '''
    Raise an exception unless all passed directories exist.
    '''

    for dirname in dirnames:
        if not is_dir(dirname):
            raise BetseDirException(
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
            raise BetseDirException(
                'Directory "{}" already exists.'.format(dirname))

# ....................{ TESTERS                            }....................
@type_check
def is_dir(dirname: str) -> bool:
    '''
    `True` only if the passed directory exists.
    '''

    return path.isdir(dirname)

# ....................{ GETTERS                            }....................
def get_current_dirname() -> str:
    '''
    Get the **current working dirname** (i.e., absolute path of the current
    working directory (CWD)) of the active Python process.

    Unless subsequently changed, this is the absolute path of the directory from
    which BETSE was initially run.
    '''

    return os.getcwd()

# ....................{ SETTERS                            }....................
@type_check
def set_current(dirname: str) -> None:
    '''
    Set the **current working directory** (CWD) of the active Python process to
    the passed directory.

    This function permanently changes the CWD for the remainder of this process.
    For a robust alternative changing the CWD for a single code block, consider
    using the :func:`current` context manager instead.

    Parameters
    -----------
    dirname : str
        Relative or absolute path of the directory to change to.
    '''

    # Log this change.
    logs.log_debug('Changing current working directory to "%s".', dirname)

    # Change to this directory.
    os.chdir(dirname)

# ....................{ CONTEXTS                           }....................
@contextmanager
@type_check
def current(dirname: str) -> GeneratorType:
    '''
    Context manager setting the **current working directory** (CWD) of the
    active Python process to the passed directory for the duration of this
    context.

    This context manager guaranteeably reverts the CWD to the prior CWD even
    when fatal exceptions are raised (e.g., due to this directory not existing).

    Parameters
    -----------
    dirname : str
        Relative or absolute path of the directory to change to.

    Returns
    -----------
    contextlib._GeneratorContextManager
        Context manager changing the CWD as described above.

    Yields
    -----------
    None
        Since this context manager yields no value, the caller's `with`
        statement must be suffixed by _no_ `as` clause.

    See Also
    -----------
    https://stackoverflow.com/a/24176022/2809027
        StackOverflow answer strongly inspiring this implementation.

    Examples
    -----------
    >>> from betse.util.paths import dirs
    >>> print('CWD: ' + dirs.get_current_dirname())
    CWD: /home/azrael
    >>> with dirs.current('/home/uriel/urial/nuriel/uryan/jeremiel'):
    ...     print('CWD: ' + dirs.get_current_dirname())
    ...     raise ValueError(
    ...         'But unknown, abstracted, brooding secret the dark power hid')
    CWD: /home/uriel/urial/nuriel/uryan/jeremiel
    ValueError: But unknown, abstracted, brooding secret the dark power hid.
    >>> print('CWD: ' + dirs.get_current_dirname())
    CWD: /home/azrael
    '''

    # Absolute path of the current CWD.
    dirname_prior = get_current_dirname()

    # Temporarily change to the passed directory. Since Python performs this
    # change only if this call raises no exceptions, this call need *NOT* be
    # embedded in the "try" block below.
    set_current(dirname)

    # Yield control to the body of the caller's "with" block.
    try:
        yield
    # Revert to the prior CWD even if that block raised an exception.
    finally:
        os.chdir(dirname_prior)

# ....................{ LISTERS                            }....................
def list_basenames(dirname: str) -> list:
    '''
    Get a list of the basenames of all files and subdirectories in the passed
    directory.
    '''

    die_unless_dir(dirname)
    return os.listdir(dirname)

# ....................{ MAKERS                             }....................
def canonicalize_and_make_unless_dir(dirname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed
    directory *and* create this directory if this directory does not already
    exist.

    This convenience function simply passes this directory to the
    func:`paths.canonicalize` and :func:`make_unless_dir` functions (in that
    order) and then returns this directory.

    See Also
    -----------
    :func:`make_unless_dir`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import paths
    dirname = paths.canonicalize(dirname)
    make_unless_dir(dirname)
    return dirname


@type_check
def make_unless_dir(dirname: str) -> None:
    '''
    Create the passed directory if this directory does *not* already exist.

    All nonexistent parents of this directory will also be recursively created,
    mimicking the action of the standard ``mkdir -p`` shell command.
    '''

    # If this directory does *NOT* already exist, create this directory. To
    # support logging, this condition is explicitly tested for. To avoid race
    # conditions (e.g., in the event this directory is created between testing
    # and creating this directory), we preserve the makedirs() keyword argument
    # "exist_ok=True" below.
    if not is_dir(dirname):
        # Log this creation.
        logs.log_debug('Creating directory: %s', dirname)

        # Create this directory if still needed.
        os.makedirs(dirname, exist_ok=True)


def make_parent_unless_dir(*pathnames: str) -> None:
    '''
    Create the parent directory of each passed path for any such directory that
    does *not* already exist.

    All nonexistent parents of each such directory will also be recursively
    created, mimicking the action of the standard ``mkdir -p`` shell command.
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
@type_check
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

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Copy us up the directory bomb.
    basename_source = paths.get_basename(dirname_source)
    copy(dirname_source, paths.join(dirname_target, basename_source))


@type_check
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

    # Log this copy.
    logs.log_debug(
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
        src=dirname_source,
        dst=dirname_target,
        symlinks=True,
    )
