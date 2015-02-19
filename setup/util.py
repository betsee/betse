#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Error-handling functions for `betse`-specific `setuptools` commands.'''

# ....................{ IMPORTS                            }....................
from distutils.errors import (
    DistutilsExecError, DistutilsFileError, DistutilsPlatformError
)
from os import path
from setuptools.cmd import Command
import os, shutil, subprocess, sys

# ....................{ EXCEPTIONS                         }....................
def die_if_os_non_posix():
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
def die_unless_file(filename: str, exception_message: str = None):
    '''
    Raise a fatal exception unless the passed non-directory file exists.
    '''
    assert isinstance(filename, str),\
        '"{}" not a string.'.format(filename)

    # If such file is not found, fail.
    if not path.isfile(filename):
        # If no such message was passed, default such message.
        if not exception_message:
             exception_message = 'File "{}" not found.'.format(filename)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise DistutilsFileError(exception_message)

def die_unless_symlink(filename: str):
    '''
    Raise a fatal exception unless the passed symbolic link exists.
    '''
    assert isinstance(filename, str),\
        '"{}" not a string.'.format(filename)

    if not path.islink(filename):
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
    subprocess.check_call(*command_words)

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
    return shutil.which(command_basename) is None

# ....................{ OUTPUTTERS                         }....................
def output_warning(*warnings) -> None:
    '''
    Print the passed warning message(s) to standard error.
    '''
    print('WARNING: ', *warnings, file = sys.stderr)

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

def entry_points(command_object: Command):
    '''
    Generator yielding a 3-tuple describing each script installed by the
    dictionary of `setuptools` options in the passed `setuptools` command
    object.

    Such dictionary of `setuptools` options *must* be a private attribute
    `_setup_options` of such command object.

    Such 3-tuple consists of each installed script's (in order):

    * Basename (e.g., `betse`).
    * Type string, guaranteed to be either:
      * `console` if such script is console-specific.
      * `gui` otherwise.
    * `:`-delimited entry point (e.g., `betse.cli.clicli:main`).
    '''
    assert isinstance(command_object, Command),\
        '"{}" not a setuptools command.'.format(command_object)

    for script_category, script_basename_to_entry_point\
        in command_object._setup_options['entry_points']:
        for script_basename, entry_point in\
            script_basename_to_entry_point.items():
            # Strip the suffix "_scripts" from such category.
            script_type = script_category[:-len('_scripts')]

            # Yield such 3-tuple.
            yield script_basename, script_type, entry_point

# --------------------( WASTELANDS                         )--------------------
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
