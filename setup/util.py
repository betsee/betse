#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Error-handling functions for `betse`-specific `setuptools` commands.'''

# ....................{ IMPORTS                            }....................
from distutils.errors import DistutilsExecError, DistutilsFileError
from os import path
from setuptools import Command
import importlib, os, platform, pkg_resources, shutil, subprocess, sys, time

# ....................{ EXCEPTIONS ~ command               }....................
def die_unless_command_succeeds(*command_words) -> None:
    '''
    Raise an exception unless running the passed command succeeds.

    For portability, such command *must* be passed as a list of shell words
    whose first element is the pathname of such command and all subsequent
    elements the command-line arguments to be passed to such command (e.g.,
    `['ls', '/']`).
    '''
    # If the first passed shell word is *NOT* pathable, raise an exception.
    die_unless_pathable(command_words[0])

    # Print the command to be run before doing so.
    print('Running "{}".'.format(' '.join(command_words)))

    # Keyword arguments to be passed to subprocess.check_call() below, which
    # accepts all arguments accepted by subprocess.Popen.__init__().
    popen_kwargs = {}

    # If the current platform is vanilla Windows, permit such command to inherit
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
    # handles will differ, in which case such behaviour prevents interaction
    # between the current shell and the vanilla Windows command to be run below.
    # In particular, all output from such command will be squelched.
    #
    # If at least one of stdin, stdout, or stderr are redirected to a blocking
    # pipe, setting "close_fds" to False can induce deadlocks under certain
    # edge-case scenarios. Since all such file handles default to None and hence
    # are *NOT* redirected in this case, "close_fds" may be safely set to False.
    #
    # On all other platforms, if "close_fds" is True then no handles *EXCEPT*
    # stdin, stdout, and stderr will be inherited by the child process. Hence,
    # this function fundamentally differs in subtle (and only slightly
    # documented ways) between vanilla Windows and all other platforms. Such
    # discrepancies appear to be harmful but probably unavoidable, given the
    # philosophical gulf between vanilla Windows and all other platforms.
    if is_os_windows_vanilla():
        popen_kwargs['close_fds'] = False

    # Run such command.
    subprocess.check_call(command_words, **popen_kwargs)

def die_unless_pathable(command_basename: str, exception_message: str = None):
    '''
    Raise an exception with the passed message if the passed **pathable** (i.e.,
    external command in the current `${PATH}`) does _not_ exist.

    If such pathable contains a directory separator, an exception is raised.
    '''
    # If such pathable is not found, raise an exception.
    if not is_pathable(command_basename):
        # If no such message was passed, default such message.
        if not exception_message:
             exception_message =\
                 'Command "{}" not found in the current ${{PATH}} or found but not an executable file.'.format(
                    command_basename)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsExecError(exception_message)

# ....................{ EXCEPTIONS ~ path                  }....................
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

        # Raise such exception.
        raise DistutilsFileError(exception_message)

def die_unless_file_or_not_found(
    pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path is either an existing non-special
    file *or* does not exist (e.g., if such path is an existing directory).
    '''
    # If such path exists and is *NOT* an existing non-special file, fail.
    if is_path(pathname) and not is_file(pathname):
        # If no such message was passed, default such message.
        if not exception_message:
            if is_dir(pathname):
                exception_message =\
                    'File "{}" already an existing directory.'.format(pathname)
            elif is_symlink(pathname):
                exception_message =\
                    'File "{}" already an existing symbolic link.'.format(
                        pathname)
            else:
                exception_message = 'Path "{}" not a file.'.format(pathname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsFileError(exception_message)

def die_unless_path(pathname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed path exists.
    '''
    # If such path is not found, fail.
    if not is_path(pathname):
        # If no such message was passed, default such message.
        if not exception_message:
            exception_message = 'Path "{}" not found.'.format(pathname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsFileError(exception_message)

def die_unless_dir(dirname: str, exception_message: str = None) -> None:
    '''
    Raise an exception unless the passed directory exists.
    '''
    # If such dir is not found, fail.
    if not is_dir(dirname):
        # If no such message was passed, default such message.
        if not exception_message:
            exception_message = 'Directory "{}" not found.'.format(dirname)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception. Since there exists no
        # DistutilsDirError(), we raise the next best thing.
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
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsFileError(exception_message)

def die_unless_symlink(filename: str) -> None:
    '''
    Raise an exception unless the passed symbolic link exists.
    '''
    assert isinstance(filename, str),\
        '"{}" not a string.'.format(filename)

    if not is_symlink(filename):
        raise DistutilsFileError(
            'Symbolic link "{}" not found.'.format(filename))

# ....................{ TESTERS ~ os                       }....................
def is_os_linux() -> bool:
    '''
    `True` if the current operating system is Linux.
    '''
    return platform.system() == 'Linux'

def is_os_posix() -> bool:
    '''
    `True` if the current operating system does _not_ comply with POSIX
    standards (e.g., as required for POSIX-style symbolic link support).

    Typically, this implies this system to _not_ be vanilla Microsoft Windows.
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
def is_path(pathname: str) -> bool:
    '''
    `True` if the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.exists(pathname)

def is_dir(pathname: str) -> bool:
    '''
    `True` if the passed directory exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.isdir(pathname)

def is_file(pathname: str) -> bool:
    '''
    `True` if the passed path is an existing non-directory file exists *after*
    following symbolic links.

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
    True if the external command with the passed basename exists.

    Specifically, return True if such basename is that of an executable file in
    the current `${PATH}`. If such basename contains a directory separator and
    is hence *not* a basename, an exception is raised.
    '''
    assert isinstance(command_basename, str),\
        '"{}" not a string.'.format(command_basename)

    # If such pathname is *NOT* a basename, such pathname erroneously contains a
    # directory separator. In such case, fail.
    if command_basename != path.basename(command_basename):
        raise DistutilsExecError(
            '"{}" contains a directory separator.'.format(command_basename))

    # Return whether such command is found.
    return shutil.which(command_basename) is not None

# ....................{ TESTERS ~ module                   }....................
def is_module(module_name: str) -> bool:
    '''
    `True` if the module with the passed fully-qualified name is importable
    under the active Python interpreter.

    If this module is a **submodule** (i.e., contains a `.` character), all
    parent modules of this module will be imported as a side effect of this
    function call. Likewise, if this module is _not_ importable via standard
    mechanisms (e.g., the OS X-specific `PyObjCTools` package), the module
    itself may also be imported as a side effect.
    '''
    # See betse.util.python.modules.is_module() for implementation details.
    assert isinstance(module_name, str),\
        '"{}" not a string.'.format(module_name)
    try:
        return importlib.util.find_spec(module_name) is not None
    except ValueError:
        try:
            importlib.import_module(module_name)
            return True
        except ImportError:
            return False
        # print('ValueError!')
        # return False

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
#FIXME: Sufficiently useful that we should probably copy this, once complete,
#to a new betse.util.io.commands.get_output() function.

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
	newlines (as under most POSIX shells), _and_ decoded via the current
        locale's preferred encoding (e.g., UTF-8).
    '''
    command_output = subprocess.check_output(args,
        # Redirect standard error to output.
        stderr = subprocess.STDOUT,

        # Decode such output via the current locale's preferred encoding.
        universal_newlines = True,
    )

    # Get such output, stripped of all trailing newlines.
    return command_output.rstrip('\n')

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
    basename *not* including such `.`) of the passed path if such path has a
    filetype or `None` otherwise.

    If such path has multiple filetypes (e.g., `odium.reigns.tar.gz`), only the
    last such filetype will be returned.
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
    print('WARNING:', *warnings, file = sys.stderr)

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

    If such target is an existing symlink, such symlink will be implicitly
    removed before being recreated.

    If such source does _not_ exist, an exception will be raised. Hence, this
    function does _not_ support creation of **dangling symbolic links** (i.e.,
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
    # If such link does *NOT* exist, fail.
    die_unless_symlink(filename)

    # Remove such link. Since symbolic links are special files, remove_file()
    # fails when passed such link and hence must be reimplemented here.
    print('Removing symbolic link "{}".'.format(filename))
    os.unlink(filename)

# ....................{ SETUPTOOLS                         }....................
def add_setup_command_classes(
    metadata: dict, setup_options: dict, *command_classes) -> None:
    '''
    Add one application-specific `setuptools` command for each passed class to
    the passed dictionary of `setuptools` options.

    For simplicity, the name of each such command will be the name of the
    corresponding class. Hence, the names of such classes are recommended to be
    short lowercase strings (e.g., `freeze`, `symlink`).
    '''
    assert isinstance(metadata, dict),\
        '"{}" not a dictionary.'.format(metadata)
    assert isinstance(setup_options, dict),\
        '"{}" not a dictionary.'.format(setup_options)

    # Add each such command class as a new command of the same name.
    for command_class in command_classes:
        assert isinstance(command_class, type),\
            '"{}" not a class.'.format(command_class)

        # Add such command.
        setup_options['cmdclass'][command_class.__name__] = command_class

        # Expose the passed dictionaries to such class by monkey-patching
        # correspoding private class fields into such classes. While passing
        # such dictionaries to instances of such class (e.g., on initialization)
        # would be ideal, distutils and hence setuptools requires commands to be
        # added as classes rather than instances. (Thus, the current approach.)
        command_class._metadata = metadata
        command_class._setup_options = setup_options

# ....................{ SETUPTOOLS ~ wrappers : generators }....................
def command_entry_points(command: Command):
    '''
    Generator yielding a 3-tuple detailing each wrapper script installed for the
    **Python distribution** (i.e., top-level package) identified by the passed
    `setuptools` command.

    See Also
    ----------
    `dist_entry_points`
        For further details on tuple contents.
    '''
    assert isinstance(command, Command),\
        '"{}" not a setuptools command.'.format(command)

    # Make a "pkg_resources"-specific distribution from the passed command. Yes,
    # this code was ripped wholesale from the run() method defined by module
    # "setuptools.command.install_scripts". Yes, we don't know how it works.
    # "Frankly, Mam, we don't give a damn."
    #
    # It should be noted that all commands have an attribute "distribution".
    # Naturally, this is a setuptools-specific distribution that has literally
    # *NOTHING* to do with the "pkg_resources"-specific distribution.
    #
    # Die, setuptools. Die!
    ei_cmd = command.get_finalized_command('egg_info')
    distribution = pkg_resources.Distribution(
        ei_cmd.egg_base,
        pkg_resources.PathMetadata(ei_cmd.egg_base, ei_cmd.egg_info),
        ei_cmd.egg_name,
        ei_cmd.egg_version,
    )

    # Defer to the generator provided by such function.
    yield from package_distribution_entry_points(distribution)

def package_distribution_entry_points(distribution: pkg_resources.Distribution):
    '''
    Generator yielding a 3-tuple detailing each wrapper script installed for the
    passed `pkg_resources`-specific distribution identifying a top-level Python
    package.

    Such 3-tuple consists of each such script's (in order):

    * Basename (e.g., `betse`). On both vanilla and Cygwin Microsoft Windows,
      such basename will be suffixed by the `.exe` filetype. On all other
      platforms, such basename will have no filetype.
    * Type string, guaranteed to be either:
      * `console` if such script is console-specific.
      * `gui` otherwise.
    * `EntryPoint` object, whose attributes specify the module to be imported
      and function to be run by such script.
    '''
    assert isinstance(distribution, pkg_resources.Distribution),\
        '"{}" not a setuptools distribution.'.format(distribution)

    # For each type of script wrapper...
    for script_type in 'console', 'gui':
        script_type_group = script_type + '_scripts'

        # for each script of such type...
        for script_basename, entry_point in\
            distribution.get_entry_map(script_type_group).items():
            # If the current platform is Windows and such script's basename has
            # no filetype, suffix such basename by ".exe".
            if is_os_windows() and get_path_filetype(script_basename) is None:
                script_basename += '.exe'

            # Yield such 3-tuple.
            yield script_basename, script_type, entry_point

# --------------------( WASTELANDS                         )--------------------
# def is_file(filename: str) -> bool:
#     '''
#     `True` if the passed non-special file exists.
#
#     `False` is returned if the passed file exists but is **special** (e.g.,
#     directory, symbolic link).
#     '''
#     assert isinstance(filename, str), '"{}" not a string.'.format(filename)
#     return path.isfile(filename)

    # Python, you win the balls.
    # return sys.path[0]
    # assert isinstance(pathname_source, str),\
    #     '"{}" not a string.'.format(pathname_source)
    # assert len(), 'Source pathname empty.'
# ....................{ EXCEPTIONS ~ os                    }....................
# def die_if_os_non_posix() -> None:
#     '''
#     Raise a fatal exception if the current operating system does `not` comply
#     with POSIX standards (e.g., as required for symbolic link manipulation).
#
#     Typically, this implies such system to be Windows.
#     '''
#     if not is_os_posix():
#         raise DistutilsPlatformError(
#             'This command requires POSIX compatibility.\n'
#             'However, the current operating system is POSIX-incompatible (e.g., Windows).'
#         )

    # return platform.system().startswith('CYGWIN_NT-')
    # return platform.system() == 'Windows'

    # The latter constraint implies shell quo this function to *not* be a general-purpose  inherently
    # If such path is *NOT* a symbolic link, fail.
    # Remove such link.
    # print('Removing symbolic link "{}".'.format(filename))
    # os.unlink(filename)

# ....................{ GETTERS                            }....................
# def get_os_type() -> bool:
#     '''
#     Get the type of the current operating system as a human-readable word.
#
#     Specifically, this is:
#
#     * `Darwin` if such system is Apple OS X.
#     * `Linux` if such system is a Linux distribution.
#     * `Windows` if such system is Microsoft Windows.
#     '''
#     return platform.system()

    # # Iterate script types.
    # for script_type in 'console', 'gui':
    #     script_type_group = script_type + '_scripts'
    #
    #     # Yield such 3-tuple for each script of such type.
    #     for script_basename, entry_point in\
    #         command.dist.get_entry_map(script_type_group).items():
    #         yield script_basename, script_type, entry_point
    # * `:`-delimited entry point (e.g., `betse.cli.clicli:main`).
#             for entry_point_spec in entry_point_specs:
#                 # Basename (e.g., "betse") and entry point (e.g.,
#                 # "betse.cli.cli:main") of such script split from such
#                 # specification on "=", stripping all leading and trailing
#                 # whitespace from such substrings after doing so.
#                 script_basename, entry_point = entry_point_spec.split('=')
#                 script_basename = script_basename.strip()
#                 entry_point = entry_point.strip()

    # for script_category, script_basename_to_entry_point\
    #     in command_object._setup_options['entry_points'].items():
    #     for script_basename, entry_point in\
    #         script_basename_to_entry_point.items():
    #         # Strip the suffix "_scripts" from such category.
    #         script_type = script_category[:-len('_scripts')]
    #
    #         # Yield such 3-tuple.
    #         yield script_basename, script_type, entry_point
    #
    # Generator yielding a 3-tuple detailing each script installed by the
    # dictionary of `setuptools` options provided by the passed `setuptools`
    # command's private attribute `_setup_options`.
#Such dictionary *must* be a private attribute `_setup_options` of such
    # command object.
    # assert isinstance(command_words, list),\
    #     '"{}" not a list.'.format(command_words)

# def die_unless_command_succeeds(
#     command_words: list,
#     *args, **kwargs) -> None:
#     '''
#     Raise an exception unless running the passed command succeeds, passing all
#     passed arguments and keyword arguments to `subprocess.call()`.
#
#     The first passed argument *must* be the command to be run. For portability,
#     such command should typically be specified as a list of shell words whose
#     first element is the pathname of such command and all subsequent elements
#     the command-line arguments to be passed such command (e.g., `['ls', '/']`).
#     If the keyword argument `shell = True` is also subsequently passed *and* a
#     Unix shell is available to the active Python3 interpreter (thus excluding
#     non-Cygwin Windows), such command may also be specified as a simple string.
#     '''
    # die_unless_pathable(command)
    # If such command fails, a CalledProcessError exception is raised.
#FUXME: Common functionality. Contemplate a utility function. To implement
#such function, we'd probably want such function to accept a list of class
#objects rather than class names, and then simply access the "__name__"
#attribute of each class object to dynamically obtain its name. Quite a bit
#simpler than the current approach, when one considers it.

# from setuptools.cmd import Command
#FUXME: Actually, this function currently only searches for command basenames in
#the current ${PATH} -- the most common usage. *shrug*

# def die_unless_pathable(command_name: str):
#     '''
#     Raise a fatal exception if the external command with the passed pathname
#     either does *not* exist *or* does but is *not* executable.
#
#     Such command will be searched for in a manner dependent on such pathname.
#     Specifically, if such pathname is:
#
#     * A basename (e.g., `ls`), the current `${PATH}` will be searched.
#     * A relative filename (e.g., `./ls`), such filename's dirname will be
#       searched relative to the current directory.
#     * An absolute filename (e.g., `/usr/bin/ls`), such filename will be
#       searched for as is.
#     '''
