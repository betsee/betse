#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level Python facilities.

This module provides functionality pertaining to the active Python interpreter
as a whole.

Caveats
----------
Word size-specific functions (e.g., :func:`is_wordsize_64`) are generally
considered poor form. Call these functions *only* where necessary.
'''

# ....................{ IMPORTS                            }....................
import platform, sys
from betse import metadata
from betse.exceptions import BetsePyException
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, MappingOrNoneTypes, SequenceTypes)

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Validate the active Python interpreter.

    This function does *not* validate this interpreter's version, as the
    top-level :mod:`betse.metadata` submodule already does so at the start of
    application startup. Instead, this function (in order):

    #. Logs a non-fatal warning if this interpreter is *not* 64-bit.
    '''

    # If this Python interpreter is 32- rather than 64-bit, log a non-fatal
    # warning. While technically feasible, running BETSE under 32-bit Python
    # interpreters imposes non-trivial constraints detrimental to sanity.
    if is_wordsize_32():
        logs.log_warning(
            '32-bit Python interpreter detected. '
            '{name} will be confined to low-precision datatypes and '
            'at most 4GB of available RAM, '
            'impeding the reliability and scalability of modelling. '
            'Consider running {name} only under a '
            '64-bit Python interpreter.'.format(name=metadata.NAME)
        )

# ....................{ TESTERS                            }....................
def is_testing() -> bool:
    '''
    ``True`` only if the active Python interpreter is running a test session
    (e.g., with the :mod:`pytest` test harness).
    '''

    return metadata._IS_TESTING

# ....................{ TESTERS ~ wordsize                 }....................
def is_wordsize_32() -> bool:
    '''
    ``True`` only if the active Python interpreter is **32-bit** (i.e., was
    compiled with a 32-bit toolchain into a 32-bit executable).
    '''

    return not is_wordsize_64()


def is_wordsize_64() -> bool:
    '''
    ``True`` only if the active Python interpreter is **64-bit** (i.e., was
    compiled with a 64-bit toolchain into a 64-bit executable).
    '''

    # Avoid circular import dependencies.
    from betse.util.type.numeric import ints

    # Return whether or not the maximum integer size supported by this Python
    # interpreter is larger than the maximum value for variables of internal
    # type `Py_ssize_t` under 32-bit Python interpreters. While somewhat obtuse,
    # this condition is well-recognized by the Python community as the optimal
    # means of portably testing this. Less preferable alternatives include:
    #
    # * "return 'PROCESSOR_ARCHITEW6432' in os.environ", which depends upon
    #   optional environment variables and hence is clearly unreliable.
    # * "return platform.architecture()[0] == '64bit'", which fails under:
    #   * OS X, returning "64bit" even when the active Python interpreter is a
    #     32-bit executable binary embedded in a so-called "universal binary."
    return sys.maxsize > ints.INT_VALUE_MAX_32_BIT

# ....................{ GETTERS                            }....................
def get_wordsize() -> int:
    '''
    Size in bits of variables of internal type `Py_ssize_t` for the active
    Python interpreter.

    This function is guaranteed to return:

    * `64`, if this is a 64-bit interpreter.
    * `32`, if this is a 32-bit interpreter.
    '''

    return 64 if is_wordsize_64() else 32


def get_version() -> str:
    '''
    Human-readable ``.``-delimited string specifying the most recent version of
    the Python language supported by the active Python interpreter (e.g.,
    ``2.7.10``, ``3.4.1``).
    '''

    return platform.python_version()

# ....................{ GETTERS ~ filename                 }....................
def get_command_line_prefix() -> list:
    '''
    List of one or more shell words unambiguously running the executable binary
    for the active Python interpreter and machine architecture.

    Since the absolute path of the executable binary for the active Python
    interpreter is insufficient to unambiguously run this binary under the
    active machine architecture, this function should *always* be called in lieu
    of :func:`get_filename` when attempting to rerun this interpreter as a
    subprocess of the current Python process. As example:

    * Under macOS, the executable binary for this interpreter may be bundled
      with one or more other executable binaries targetting different machine
      architectures (e.g., 32-bit, 64-bit) in a single so-called "universal
      binary." Distinguishing between these bundled binaries requires passing
      this interpreter to a prefixing macOS-specific command, ``arch``.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # List of such shell words.
    command_line = None

    # If this is OS X, this interpreter is only unambiguously runnable via the
    # OS X-specific "arch" command.
    if oses.is_macos():
        # Run the "arch" command.
        command_line = ['arch']

        # Instruct this command to run the architecture-specific binary in
        # Python's universal binary corresponding to the current architecture.
        if is_wordsize_64():
            command_line.append('-i386')
        else:
            command_line.append('-x86_64')

        # Instruct this command, lastly, to run this interpreter.
        command_line.append(get_filename())
    # Else, this interpreter is unambiguously runnable as is.
    else:
        command_line = [get_filename()]

    # Return this list.
    return command_line


def get_filename() -> str:
    '''
    Absolute filename of the executable binary for the active Python process.
    '''

    # Absolute filename of this executable if this path is retrievable by Python
    # *OR* either "None" or the empty string otherwise.
    py_filename = sys.executable

    # If this filename is *NOT* retrievable, raise an exception.
    if not py_filename:
        raise BetsePyException(
            'Absolute filename of Python interpreter not retrievable.')

    # Return this filename.
    return py_filename


def get_shebang() -> str:
    '''
    POSIX-compliant **shebang** (i.e., machine-readable ``#!``-prefixed first
    line of interpretable scripts) unambiguously running the executable binary
    for the active Python interpreter and machine architecture.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # List of one or more shell words unambiguously running this binary.
    py_command_args = get_command_line_prefix()

    # Space-delimited string listing these words.
    py_command = strs.join_on_space(py_command_args)

    # Return this string prefixed by the shebang identifier.
    return '#!' + py_command

# ....................{ GETTERS ~ metadata                 }....................
def get_metadata() -> 'OrderedArgsDict':
    '''
    Ordered dictionary synopsizing the active Python interpreter.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyfreeze
    from betse.util.type.mapping.mapcls import OrderedArgsDict

    # Return this dictionary.
    return OrderedArgsDict(
        'version', get_version(),
        'wordsize', get_wordsize(),
        'is frozen', pyfreeze.is_frozen(),
    )

# ....................{ ADDERS                             }....................
@type_check
def add_import_dirname(dirname: str) -> None:
    '''
    Register all files and subdirectories of the directory with the passed path
    to be importable modules and packages (respectively) for the remainder of
    the current Python process if this directory has not already been registered
    *or* noop otherwise.

    Specifically, this function appends this dirname to the current
    :data:`sys.path` listing (in order) the dirnames of all directories to be
    iteratively searched for any module or package on first importing that
    module or package. To comply with Python standards in which the first item
    of this list is either the dirname of the directory containing the script
    from which this process was invoked *or* the empty string (signifying the
    current directory), this list is appended to rather than prepended to.

    Parameters
    ----------
    dirname : str
        Absolute or relative path of the directory to be registered.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Log this addition.
    logs.log_debug('Registering import directory: %s', dirname)

    # If this directory does *NOT* exist or is unreadable, raise an exception.
    dirs.die_unless_dir(dirname)

    # If the current PYTHONPATH already contains this directory...
    if dirname in sys.path:
        # Log this edge case.
        logs.log_debug('Ignoring already registered import directory.')

        # Noop.
        return
    # Else, the current PYTHONPATH does *NOT* already contain this directory.

    # Append this directory to the current PYTHONPATH.
    sys.path.append(dirname)

# ....................{ RUNNERS                            }....................
@type_check
def rerun_or_die(
    command_args: SequenceTypes,
    popen_kwargs: MappingOrNoneTypes = None) -> None:
    '''
    Rerun the active Python interpreter as a subprocess of the current Python
    process, raising an exception on subprocess failure.

    Parameters
    ----------
    command_args : SequenceTypes
        List of zero or more arguments to pass to this interpreter.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    See Also
    ----------
    :func:`betse.util.path.command.cmdrun.run_or_die`
        Low-level commentary on subprocess execution.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import cmdrun

    # List of one or more shell words comprising this command.
    command_words = get_command_line_prefix() + command_args

    # Rerun this interpreter.
    return cmdrun.run_or_die(
        command_words=command_words, popen_kwargs=popen_kwargs)
