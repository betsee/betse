#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
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

* Assumes BETSE to be available. While this is certainly the case when this file
  resides in the BETSE codebase, this is *not* necessarily the case when this
  file is copied into and hence resides in the codebases of other projects
  (e.g., BETSEE). In these projects, BETSE is merely yet another dependency that
  is typically unavailable at installation time.
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

# ....................{ IMPORTS                            }....................
import importlib, os, platform, shutil, subprocess, sys, time
from distutils.errors import (
    DistutilsExecError, DistutilsFileError, DistutilsModuleError)
from os import path
from pkg_resources import Distribution, PathMetadata
from setuptools import Command

# ....................{ EXCEPTIONS ~ command               }....................
def die_unless_command_succeeds(*command_words) -> None:
    '''
    Raise an exception unless running the passed command succeeds.

    For portability, this command *must* be passed as a list of shell words
    whose first element is the pathname of such command and all subsequent
    elements the command-line arguments to be passed to such command (e.g.,
    ``['ls', '/']``).
    '''

    # If the first passed shell word is *NOT* pathable, raise an exception.
    die_unless_pathable(command_words[0])

    # Print the command to be run before doing so.
    print('Running "{}".'.format(' '.join(command_words)))

    # Keyword arguments to be passed to subprocess.check_call() below, which
    # accepts all arguments accepted by subprocess.Popen.__init__().
    popen_kwargs = {}

    # If the current platform is vanilla Windows, permit this command to inherit
    # all file handles (including stdin, stdout, and stderr) from the current
    # process. By default, subprocess.Popen documentation insists that:
    #
    #     "On Windows, if close_fds is true then no handles will be inherited by
    #      the child process."
    #
    # The child process will then open new file handles for stdin, stdout, and
    # stderr. If the current terminal is a Windows Console, the underlying
    # terminal devices and hence file handles will remain the same, in which
    # case this is *NOT* an issue. If the current terminal is Cygwin-based
    # (e.g.,, MinTTY), however, the underlying terminal devices and hence file
    # handles will differ, in which case this behaviour prevents interaction
    # between the current shell and the vanilla Windows command to be run below.
    # In particular, all output from this command will be squelched.
    #
    # If at least one of stdin, stdout, or stderr are redirected to a blocking
    # pipe, setting "close_fds" to False can induce deadlocks under certain
    # edge-case scenarios. Since all such file handles default to None and hence
    # are *NOT* redirected in this case, "close_fds" may be safely set to False.
    #
    # On all other platforms, if "close_fds" is True, no file handles *EXCEPT*
    # stdin, stdout, and stderr will be inherited by the child process. Hence,
    # this function fundamentally differs in subtle (and only slightly
    # documented ways) between vanilla Windows and all other platforms. These
    # discrepancies appear to be harmful but probably unavoidable, given the
    # philosophical gulf between vanilla Windows and all other platforms.
    if is_os_windows_vanilla():
        popen_kwargs['close_fds'] = False

    # Run this command.
    subprocess.check_call(command_words, **popen_kwargs)


def die_unless_pathable(command_basename: str, exception_message: str = None):
    '''
    Raise an exception with the passed message if the passed **pathable** (i.e.,
    external command in the current `${PATH}`) does _not_ exist.
    '''

    # If this pathable is not found, raise an exception.
    if not is_pathable(command_basename):
        # If no message was passed, default this message.
        if not exception_message:
             exception_message = (
                 'Command "{}" not found in the current ${{PATH}} or '
                 'found but not an executable file.'.format(
                    command_basename))
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsExecError(exception_message)

# ....................{ EXCEPTIONS ~ path                  }....................
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
                    'File "{}" already an existing directory.'.format(pathname))
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

# ....................{ EXCEPTIONS ~ python                }....................
def die_unless_module(module_name: str, exception_message: str = None):
    '''
    Raise an exception with the passed message if the module with the passed
    fully-qualified name (e.g., `astarte.ashtoreth.ishtar`) is unimportable.
    '''

    # If this module is unimportable, raise an exception.
    if not is_module(module_name):
        # If no message was passed, default this message.
        if not exception_message:
             exception_message = (
                 'Module "{}" not installed or not importable under '
                 'the current Python interpreter.'.format(module_name))
        assert isinstance(exception_message, str), (
            '"{}" not a string.'.format(exception_message))

        # Raise this exception.
        raise DistutilsModuleError(exception_message)

# ....................{ TESTERS ~ os                       }....................
def is_os_linux() -> bool:
    '''
    `True` if the current operating system is Linux.
    '''
    return platform.system() == 'Linux'


def is_os_posix() -> bool:
    '''
    `True` if the current operating system complies with POSIX standards (e.g.,
    as required for POSIX-compliant symbolic link support).

    Typically, this implies this system to _not_ be vanilla Microsoft Windows
    (i.e., to be either a Cygwin-enabled Windows terminal _or_ a genuine
    POSIX-compliant system).
    '''
    return os.name == 'posix'
    # return False


def is_os_os_x() -> bool:
    '''
    `True` if the current operating system is Apple OS X.
    '''
    return platform.system() == 'Darwin'

# ....................{ TESTERS ~ os : windows             }....................
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

# ....................{ TESTERS ~ path                     }....................
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


def is_pathable(command_basename: str) -> bool:
    '''
    `True` only if the external command with the passed basename exists (i.e.,
    is an executable file in the current `${PATH}`).
    '''

    # Sanitize this command basename.
    command_basename = sanitize_command_basename(command_basename)

    # Return whether this command is found or not.
    return shutil.which(command_basename) is not None

# ....................{ TESTERS ~ module                   }....................
def is_module(module_name: str) -> bool:
    '''
    `True` only if the module with the passed fully-qualified name is importable
    under the active Python interpreter.

    If this module is a **submodule** (i.e., contains a `.` character), all
    parent modules of this module will be imported as a side effect of this
    function call. Likewise, if this module is _not_ importable via standard
    mechanisms (e.g., the OS X-specific `PyObjCTools` package), the module
    itself may also be imported as a side effect.
    '''

    # See betse.util.python.modules.is_module() for implementation details.
    assert isinstance(module_name, str), (
        '"{}" not a string.'.format(module_name))
    assert len(module_name), 'Module name empty.'
    try:
        return importlib.util.find_spec(module_name) is not None
    except ValueError:
        try:
            importlib.import_module(module_name)
            return True
        except ImportError:
            return False

# ....................{ GETTERS                            }....................
def get_project_dirname():
    '''
    Get the absolute path of the directory containing the currently run
    `setup.py` script.
    '''
    # While such path is also typically  available as the first entry of the
    # "sys.path" list, you know what they say about assumptions.
    return get_path_dirname(get_path_dirname(__file__))

# ....................{ GETTERS ~ io                       }....................
def get_command_output(*args) -> str:
    '''
    Get all standard output and error captured by running the external shell
    command signified by the passed list if this command succeeds or raise an
    exception detailing this command's failure otherwise.

    Parameters
    ----------
    *args : list
        List of shell words comprising this command. The first item of this list
        should be the pathname for this command; all remaining items should be
        the arguments to pass this command.

    Returns
    ----------
    str
        All standard output and error captured by running this command,
        interleaved together in output order, stripped of all trailing
        newlines (as under most POSIX shells), *and* decoded via the current
        locale's preferred encoding (e.g., UTF-8).
    '''

    command_output = subprocess.check_output(
        args,

        # Redirect standard error to output.
        stderr=subprocess.STDOUT,

        # Decode this output using the current locale's preferred encoding.
        universal_newlines=True,
    )

    # Get this output, stripped of all trailing newlines.
    return command_output.rstrip('\n')

# ....................{ GETTERS ~ io : file                }....................
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

# ....................{ GETTERS ~ path                     }....................
def get_path_canonicalized(pathname: str) -> str:
    '''
    Get the **canonical form** (i.e., unique absolute path) of the passed path.

    Specifically (in order):

    * Perform **tilde expansion,** replacing a `~` character prefixing such path
      by the absolute path of the current user's home directory.
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
    Get the **dirname** (i.e., parent directory) of the passed path if such path
    has a dirname or raise an exception otherwise.
    '''
    # Get such dirname. Since the path.dirname() function returns the empty
    # string for paths containing no directory separators and hence having no
    # dirnames, assert such return value to be non-empty.
    dirname = path.dirname(pathname)
    assert len(dirname), 'Pathname "{}" dirname empty.'.format(pathname)
    return dirname

# ....................{ GETTERS ~ path : filetype          }....................
def get_path_filetype(pathname: str) -> str:
    '''
    Get the **last filetype** (i.e., last `.`-prefixed substring of the
    basename *not* including such `.`) of the passed path if this path has a
    filetype _or_ `None` otherwise.

    If this path contains multiple filetypes (e.g., `odium.reigns.tar.gz`), this
    function returns only the last filetype.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'

    # Such filetype. (Yes, splitext() is exceedingly poorly named.)
    filetype = path.splitext(pathname)[1]

    # Get such filetype, stripping the prefixing "." from the string returned by
    # the prior call if such path has a filetype or returning None otherwise.
    return filetype[1:] if filetype else None


def get_path_sans_filetype(pathname: str) -> str:
    '''
    Get the passed path without last filetype (including prefixing `.`) if such
    path has a filetype *or* as is otherwise.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    assert len(pathname), 'Pathname empty.'
    return path.splitext(pathname)[0]

# ....................{ SANITIZERS                         }....................
def sanitize_command_basename(command_basename: str) -> str:
    '''
    Convert the passed platform-agnostic command basename (e.g., `pytest`) into
    a platform-specific command basename (e.g., `pytest.exe`).

    If the passed basename contains a directory separator and hence is _not_ a
    basename, an exception is raised. Else, under:

    * Windows, the passed basename is appended by `.exe`. To avoid confusion
      with non-Windows executables in the current `${PATH}` when running under
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

# ....................{ IMPORTERS                          }....................
def import_module(module_name: str, exception_message: str = None) -> type(sys):
    '''
    Dynamically import and return the module, package, or C extension with the
    passed fully-qualified name.

    If this module is unimportable, an exception with the passed message is
    raised.
    '''

    # If this module is unimportable, raise an exception.
    die_unless_module(module_name, exception_message)

    # Else, import and return this module.
    return importlib.import_module(module_name)

# ....................{ QUITTERS                           }....................
def exit_with_status(exit_status: int) -> None:
    '''
    Terminate the current Python process with the passed 0-based exit status.
    '''
    sys.exit(exit_status)

# ....................{ OUTPUTTERS                         }....................
def output_sans_newline(*strings) -> None:
    '''
    Print the passed strings to standard output *not* suffixed by a newline.

    By default, printed strings are suffixed by a newline.
    '''
    print(*strings, end = '')


def output_warning(*warnings) -> None:
    '''
    Print the passed warning messages to standard error.
    '''
    print('WARNING: ', *warnings, file=sys.stderr)

# ....................{ QUOTERS                            }....................
def shell_quote(text: str) -> str:
    '''
    Shell-quote the passed string.

    If the current operating system is:

    * *Not* Windows (e.g., Linux, OS X), the returned string is guaranteed to be
      suitable for passing as an arbitrary positional argument to external
      commands.
    * Windows, the returned string is suitable for passing *only* to external
      commands parsing arguments according in the same manner as the Microsoft C
      runtime. Whereas *all* applications running under POSIX-compliant systems
      are required to parse arguments in the same manner (e.g., according to
      Bourne shell lexing), no such standard applies to applications running
      under Windows. For this reason, shell quoting is inherently unreliable
      under Windows.
    '''
    assert isinstance(text, str), '"{}" not a string.'.format(text)

    # If the current OS is vanilla Windows, do *NOT* perform POSIX-compatible
    # quoting. Vanilla Windows is POSIX-incompatible and hence does *NOT* parse
    # command-line arguments according to POSIX standards. In particular,
    # vanilla Windows does *NOT* treat single-quoted arguments as single
    # arguments but rather as multiple shell words delimited by the raw literal
    # `'`. This is circumventable by calling an officially undocumented
    # Windows-specific Python function. (Awesome.)
    if is_os_windows_vanilla():
        return subprocess.list2cmdline([text])
    # Else, perform POSIX-compatible quoting.
    else:
        import shlex
        return shlex.quote(text)

# ....................{ MAKERS                             }....................
def make_dir_unless_found(dirname: str) -> None:
    '''
    Create the passed directory if such directory does *not* already exist.

    All nonexistent parents of such directory will also be recursively created,
    mimicking the action of the conventional shell command `mkdir -p`.
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
        print('Creating directory "{}".'.format(dirname))

        # Create such directory if still needed.
        os.makedirs(dirname, exist_ok = True)


def make_symlink(pathname_source: str, filename_target: str) -> None:
    '''
    Symbolically link the passed source path to the passed target symlink.

    If this target is an existing symlink, this symlink will be implicitly
    removed before being recreated.

    If this source does *not* exist, an exception will be raised. Hence, this
    function does *not* support creation of **dangling symbolic links** (i.e.,
    links to non-existent paths).
    '''

    # If such source path does *NOT* exist, raise an exception.
    die_unless_path(pathname_source)

    # If such link currently exists, remove such link.
    if is_symlink(filename_target):
        remove_symlink(filename_target)

    # (Re)create such link.
    print('Symbolically linking "{}" to "{}".'.format(
        pathname_source, filename_target))
    os.symlink(pathname_source, filename_target)

# ....................{ MOVERS                             }....................
def move_file(filename_source: str, filename_target: str) -> None:
    '''
    Move the passed source to the passed target file.
    '''

    # If such file does *NOT* exist, raise an exception.
    die_unless_file(filename_source)

    # Move such file.
    print('Moving file "{}" to "{}".'.format(filename_source, filename_target))
    shutil.move(filename_source, filename_target)

# ....................{ REMOVERS                           }....................
def remove_path(pathname: str) -> None:
    '''
    Recursively remove the passed directory in a safe manner (e.g., *not*
    following symbolic links outside such directory).

    This is an inherently dangerous operation and hence delayed for several
    seconds, allowing sufficiently aware users to jam the panic button.
    '''

    # If such path does *NOT* exist, fail.
    die_unless_path(pathname)

    # If such path is a directory, remove such directory.
    if is_dir(pathname):
        remove_dir(pathname)
    # Else, remove such file.
    else:
        remove_file(pathname)


def remove_dir(dirname: str) -> None:
    '''
    Recursively remove the passed directory in a safe manner (e.g., *not*
    following symbolic links outside such directory).

    This is an inherently dangerous operation and hence delayed for several
    seconds, allowing sufficiently aware users to jam the panic button.
    '''

    # If such directory does *NOT* exist, fail.
    die_unless_dir(dirname)

    # For safety, wait several seconds to do so. (Read: panic button.)
    sleep_seconds = 4
    print('Removing directory "{}" in {} seconds...'.format(
        dirname, sleep_seconds))
    time.sleep(sleep_seconds)

    # Remove such directory.
    shutil.rmtree(dirname)
    print('Removed directory "{}".'.format(dirname))


def remove_file(filename: str) -> None:
    '''
    Remove the passed non-special file.
    '''

    # If such file does *NOT* exist, fail.
    die_unless_file(filename)

    # Remove such file.
    print('Removing file "{}".'.format(filename))
    os.unlink(filename)


def remove_symlink(filename: str) -> None:
    '''
    Remove the passed symbolic link.
    '''

    # If this link does *NOT* exist, fail.
    die_unless_symlink(filename)

    # Remove this link. Since symbolic links are special files, remove_file()
    # fails when passed such link and hence must be reimplemented here.
    print('Removing symbolic link "{}".'.format(filename))
    os.unlink(filename)

# ....................{ SETUPTOOLS                         }....................
def add_setup_command_classes(
    metadata: dict, setup_options: dict, *command_classes) -> None:
    '''
    Add one application-specific :mod:`setuptools` command for each passed class
    to the passed dictionary of :mod:`setuptools` options.

    For simplicity, the name of each such command will be the name of the
    corresponding class. Hence, the names of such classes are recommended to be
    short lowercase strings (e.g., ``freeze``, ``symlink``).
    '''
    assert isinstance(metadata, dict), '"{}" not a dictionary.'.format(metadata)
    assert isinstance(setup_options, dict), (
        '"{}" not a dictionary.'.format(setup_options))

    # Add each such command class as a new command of the same name.
    for command_class in command_classes:
        assert isinstance(command_class, type), (
            '"{}" not a class.'.format(command_class))

        # Add this command.
        setup_options['cmdclass'][command_class.__name__] = command_class

        # Expose the passed dictionaries to this class by monkey-patching
        # application-specific private class variables into these classes. While
        # passing these dictionaries to instances of this class (e.g., on
        # instantiation) would be ideal, distutils and hence setuptools
        # requires commands to be registered as classes rather than instances.
        command_class._metadata = metadata
        command_class._setup_options = setup_options

# ....................{ SETUPTOOLS ~ wrappers : generators }....................
def command_entry_points(command: Command) -> 'GeneratorType':
    '''
    Generator yielding a 3-tuple detailing each wrapper script installed for the
    **Python distribution** (i.e., top-level package) identified by the passed
    `setuptools` command.

    See Also
    ----------
    `package_distribution_entry_points()`
        Further details on tuple contents.
    '''
    assert isinstance(command, Command), (
        '"{}" not a setuptools command.'.format(command))

    # Make a "pkg_resources"-specific distribution from the passed command. Yes,
    # this code was ripped wholesale from the run() method defined by module
    # "setuptools.command.install_scripts". Yes, we don't know how it works.
    # "Frankly, Mam, we don't give a damn."
    #
    # It should be noted that all commands have an attribute "distribution".
    # Naturally, this is a setuptools-specific distribution that has literally
    # *NOTHING* to do with pkg_resources-style distributions.
    #
    # Die, setuptools. Die!
    ei_cmd = command.get_finalized_command('egg_info')
    distribution = Distribution(
        ei_cmd.egg_base,
        PathMetadata(ei_cmd.egg_base, ei_cmd.egg_info),
        ei_cmd.egg_name,
        ei_cmd.egg_version,
    )

    # Defer to the generator provided by such function.
    yield from package_distribution_entry_points(distribution)


def package_distribution_entry_points(
    distribution: '(Distribution, VersionlessRequirement)') -> None:
    '''
    Generator yielding a 3-tuple describing each wrapper script installed for
    the passed `pkg_resources`-specific distribution identifying a unique
    top-level Python package.

    Yields
    ----------
    For each such script, this generator yields a 3-tuple
    `(script_basename, ui_type, entry_point)` such that:

    * `script_basename` is this script's basename (e.g., `betse`). To simplify
      integration with the downstream setuptools API (e.g., the
      `setuptools.command.easy_install.ScriptWriter.get_script_args()` method),
      this basename is typically _not_ suffixed by a platform-specific filetype
      (e.g., `.exe` under vanilla or Cygwin Microsoft Windows).
    * `ui_type` is this script's interface type string, guaranteed to be either:
      * `console` if this script is console-specific.
      * `gui` otherwise.
    * `entry_point` is this script's `EntryPoint` object, whose attributes
      specify the module to be imported and function to be run by this script.

    Parameters
    ----------
    distribution : Distribution, VersionlessRequirement
        Distribution object identifying the top-level Python package to yield
        entry points for. Specifically, either:
        * A `Distribution` object supplied by the `install` or `symlink`
          subcommands.
        * A `VersionlessRequirement` object supplied by the `develop`
          subcommand. As the classname suggests, this object wraps the
          corresponding `Distribution` object by stripping versioning from this
          distribution's name (e.g., reducing `foo==1.0` to merely `foo`).
    '''
    # Do *NOT* bother attempting to assert the passed distribution to be an
    # instance of either the "Distribution" or "VersionlessRequirement" classes.
    # While the former is guaranteed to exist under all setuptools version, the
    # latter is *NOT*. Instead, only assert this distribution to be non-None.
    assert distribution is not None, 'Setuptools distribution expected.'

    # For each type of script wrapper...
    for script_type in 'console', 'gui':
        script_type_group = script_type + '_scripts'

        # For each script of this type...
        for script_basename, entry_point in (
            distribution.get_entry_map(script_type_group).items()):
            # Yield this 3-tuple. To simplify integration with the downstream
            # setuptools API, do *NOT* sanitize_snakecase this script's basename by
            # calling sanitize_command_basename(). Since that API already
            # implicitly suffixes this basename by ".exe", doing so here would
            # erroneously result in this basename being suffixed by ".exe.exe".
            yield script_basename, script_type, entry_point
