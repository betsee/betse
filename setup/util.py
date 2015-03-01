#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Error-handling functions for `betse`-specific `setuptools` commands.'''

# ....................{ IMPORTS                            }....................
from distutils.errors import (
    DistutilsExecError, DistutilsFileError, DistutilsPlatformError)
from os import path
from setuptools import Command
import os, pkg_resources, shutil, subprocess, sys, time

# ....................{ EXCEPTIONS                         }....................
def die_if_os_non_posix() -> None:
    '''
    Raise a fatal exception if the current operating system does `not` comply
    with POSIX standards (e.g., as required for symbolic link manipulation).

    Typically, this implies such system to be Windows.
    '''
    if os.name != 'posix':
        raise DistutilsPlatformError(
            'This command requires POSIX compliance. Distressingly, the current '
            'operating system is POSIX-noncompliant (e.g., Windows).'
        )

# ....................{ EXCEPTIONS ~ path                  }....................
def die_unless_dir_or_not_found(
    pathname: str, exception_message: str = None) -> None:
    '''
    Raise a fatal exception unless the passed path is either an existing
    directory *or* does not exist (i.e., if such path is an existing non-
    directory).
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
    Raise a fatal exception unless the passed path is either an existing non-
    special file *or* does not exist (e.g., if such path is an existing
    directory).
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

def die_unless_dir(dirname: str, exception_message: str = None) -> None:
    '''
    Raise a fatal exception unless the passed directory exists.
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
    Raise a fatal exception unless the passed non-special file exists.
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
    Raise a fatal exception unless the passed symbolic link exists.
    '''
    assert isinstance(filename, str),\
        '"{}" not a string.'.format(filename)

    if not is_symlink(filename):
        raise DistutilsFileError(
            'Symbolic link "{}" not found.'.format(filename))

# ....................{ EXCEPTIONS ~ command               }....................
def die_unless_command_succeeds(*command_words) -> None:
    '''
    Raise an exception unless running the passed command succeeds.

    For portability, such command *must* be passed as a list of shell words
    whose first element is the pathname of such command and all subsequent
    elements the command-line arguments to be passed to such command (e.g.,
    `['ls', '/']`).
    '''
    # Die unless the first passed shell word is an existing command.
    die_unless_command(command_words[0])

    # Print the command to be run before doing so.
    print('Running "{}".'.format(' '.join(command_words)))

    # Run such command.
    subprocess.check_call(command_words)

def die_unless_command(command_basename: str, exception_message: str = None):
    '''
    Raise a fatal exception with the passed message if the external command with
    the passed basename does *not* exist.

    Specifically, raise such exception if such basename is not that of an
    executable file in the current `${PATH}`. If such basename contains a
    directory separator, an exception is also raised.
    '''
    # If such command is not found, fail.
    if not is_command(command_basename):
        # If no such message was passed, default such message.
        if not exception_message:
             exception_message =\
                 'Command "{}" not found in the current PATH or found but not an executable file.'.format(
                    command_basename)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsExecError(exception_message)

# ....................{ TESTERS                            }....................
def is_path(pathname: str) -> bool:
    '''
    True if the passed path exists.
    '''
    assert isinstance(pathname, str), '"{}" not a string.'.format(pathname)
    return path.exists(pathname)

def is_dir(dirname: str) -> bool:
    '''
    True if the passed directory exists.
    '''
    assert isinstance(dirname, str), '"{}" not a string.'.format(dirname)
    return path.isdir(dirname)

def is_file(filename: str) -> bool:
    '''
    True if the passed non-special file exists.

    This function returns False if such file exists but is **special** (e.g.,
    directory, device node, symbolic link).
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    return path.isfile(filename)

def is_symlink(filename: str) -> bool:
    '''
    True if the passed symbolic link exists.

    Caveats
    ----------
    This function returns False if the passed symbolic link exists but the
    current user has insufficient privelages to follow such link. This may
    constitute a bug in the underlying `path.islink()` function.
    '''
    assert isinstance(filename, str), '"{}" not a string.'.format(filename)
    return path.islink(filename)

def is_command(command_basename: str) -> bool:
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
    print('WARNING: ', *warnings, file = sys.stderr)

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

# ....................{ MOVERS                             }....................
def move_file(filename_source: str, filename_target: str) -> None:
    '''
    Move the passed source to the passed target file.
    '''
    # If such file does *NOT* exist, fail.
    die_unless_file(filename_source)

    # Move such file.
    print('Moving file "{}" to "{}".'.format(filename_source, filename_target))
    shutil.move(filename_source, filename_target)

# ....................{ REMOVERS                           }....................
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
    sleep_seconds = 8
    print('Removing "{}" in {} seconds...'.format(dirname, sleep_seconds))
    time.sleep(sleep_seconds)

    # Remove such directory.
    print('Removing directory "{}".'.format(dirname))
    shutil.rmtree(dirname)

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
def add_setup_command_classes(setup_options: dict, *command_classes) -> None:
    '''
    Add one application-specific `setuptools` command for each passed class to
    the passed dictionary of `setuptools` options.

    For simplicity, the name of each such command will be the name of the
    corresponding class. Hence, the names of such classes are recommended to be
    short lowercase strings (e.g., `freeze`, `symlink`).
    '''
    assert isinstance(setup_options, dict),\
        '"{}" not a dictionary.'.format(setup_options)

    # Add each such command class as a new command of the same name.
    for command_class in command_classes:
        assert isinstance(command_class, type),\
            '"{}" not a class.'.format(command_class)

        # Add such command.
        setup_options['cmdclass'][command_class.__name__] = command_class

        # Expose the passed dictionary of "setuptools" options to such class by
        # adding a new private class field "_setup_options" to such class. While
        # passing such dictionary to instances of such class would be ideal,
        # distutils and hence setuptools requires commands to be added as
        # classes rather than instances. Thus, the current approach.
        command_class._setup_options = setup_options

# ....................{ SETUPTOOLS ~ entry points          }....................
def command_entry_points(command: Command):
    '''
    Generator yielding a 3-tuple detailing each wrapper script installed for the
    *Python distribution* (i.e., top-level package) identified by the passed
    `setuptools` command.

    See Also
    ----------
    dist_entry_points
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
    # *NOTHING* to do with "pkg_resources"-specific distribution.
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

    * Basename (e.g., `betse`).
    * Type string, guaranteed to be either:
      * `console` if such script is console-specific.
      * `gui` otherwise.
    * `EntryPoint` object, whose attributes specify the module to be imported
      and function to be run by such script.
    '''
    assert isinstance(distribution, pkg_resources.Distribution),\
        '"{}" not a setuptools distribution.'.format(distribution)

    # Iterate script types.
    for script_type in 'console', 'gui':
        script_type_group = script_type + '_scripts'

        # Yield such 3-tuple for each script of such type.
        for script_basename, entry_point in\
            distribution.get_entry_map(script_type_group).items():
            yield script_basename, script_type, entry_point

# --------------------( WASTELANDS                         )--------------------
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
    # die_unless_command(command)
    # If such command fails, a CalledProcessError exception is raised.
#FUXME: Common functionality. Contemplate a utility function. To implement
#such function, we'd probably want such function to accept a list of class
#objects rather than class names, and then simply access the "__name__"
#attribute of each class object to dynamically obtain its name. Quite a bit
#simpler than the current approach, when one considers it.

# from setuptools.cmd import Command
#FUXME: Actually, this function currently only searches for command basenames in
#the current ${PATH} -- the most common usage. *shrug*

# def die_unless_command(command_name: str):
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
