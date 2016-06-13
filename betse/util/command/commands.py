#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
**Command** (i.e., external executable file) facilities.

Command Words Arguments
----------
Most **runners** (i.e., functions running commands) defined by this module
accept a mandatory `command_words` argument, a list of one or more shell words
comprising this command whose

* Mandatory first element is either:
    * This command's absolute or relative path.
    * This command's basename, in which case the first command with that
      basename in the current `${PATH}` environment variable will be run. If no
      such command is found, an exception is raised.
* Optional subsequent elements are this command's arguments (in order).

`Popen()` Keyword Arguments
----------
Most **runners** (i.e., functions running commands) defined by this module
accept optional keyword arguments accepted by the `subprocess.Popen.__init__()`
constructor. Frequently passed options of interest include:

* `cwd`, the absolute path of the current working directory (CWD) from which
  this command is to be run. Defaults to the current CWD.
* `timeout`, the maximum number of milliseconds this command is to be run for.
  Commands with execution time exceeding this timeout will be mercifully killed.
  Defaults to `None`, in which case this command will run indefinitely.
'''

# ....................{ IMPORTS                            }....................
import shutil, subprocess, sys
from betse.exceptions import BetseExceptionPath
from betse.metadata import SCRIPT_NAME_CLI
from betse.util.io.log import logs
from betse.util.type import types
from subprocess import CalledProcessError, TimeoutExpired

# ....................{ GLOBALS                            }....................
_CURRENT_BASENAME = None
'''
Basename of the command originating the active Python interpreter.

This global caches the return value of the frequently called
`get_current_basename()` function, for efficiency.
'''

# ....................{ TESTERS                            }....................
def is_pathable(command_basename: str) -> bool:
    '''
    `True` only if the external command with the passed basename exists.

    This function returns `True` only if this basename is that of an executable
    file in the current `${PATH}`. If this basename contains a directory
    separator and is hence _not_ a basename, an exception is raised.
    '''
    assert types.is_str_nonempty(command_basename), (
        types.assert_not_str_nonempty(command_basename, 'Command name'))

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # If this string is *NOT* a pure basename, fail.
    if paths.is_basename(command_basename):
        raise BetseExceptionPath(
            'Command "{}" contains directory separators.'.format(
                command_basename))

    # Return whether this command exists or not.
    return shutil.which(command_basename) is not None

# ....................{ GETTERS ~ path                     }....................
def get_current_basename() -> str:
    '''
    Get the basename of the command originating the active Python interpreter.

    If this interpreter is interpreting a block of arbitrary runtime code passed
    to this interpreter on the command line via Python's `-c` option (e.g., due
    to being called by a distributed `pytest-xdist` test), this function
    unconditionally returns the basename of the BETSE CLI (e.g., `betse`) rather
    than `-c`. `-c` is a rather unexpected and non-human-readable basename!
    '''

    # If this function has already been called, return the string cached by the
    # first call to this function.
    global _CURRENT_BASENAME
    if _CURRENT_BASENAME is not None:
        return _CURRENT_BASENAME

    # Avoid circular import dependencies.
    from betse.util.path import paths

    # Raw absolute or relative path of the current command.
    _CURRENT_BASENAME = sys.argv[0]

    # If this is the non-human-readable "-c" Python interpreter option,
    # substitute this with the human-readable basename of the BETSE CLI.
    if _CURRENT_BASENAME == '-c':
        _CURRENT_BASENAME = SCRIPT_NAME_CLI
    # Else, reduce this absolute or relative path to a basename.
    else:
        _CURRENT_BASENAME = paths.get_basename(_CURRENT_BASENAME)

    # Return this basename.
    return _CURRENT_BASENAME

# ....................{ RUNNERS                            }....................
def run(command_words: list, **popen_kwargs) -> None:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure.

    This exception contains the exit status of this subprocess.

    Parameters
    ----------
    command_words : list
        List of one or more shell words comprising this command.
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Avoid circular import dependencies.
    from betse.util.command import exits

    # Run this command, raising an exception on command failure. For
    # reusability, reimplement the subprocess.check_call() function here rather
    # than explicitly call this function. The latter approach would require
    # duplicating logic between this and the run_nonfatal() function.
    exit_status = run_nonfatal(command_words, popen_kwargs)
    if exits.is_failure(exit_status):
        raise CalledProcessError(exit_status, command_words)


def run_nonfatal(command_words: list, **popen_kwargs) -> int:
    '''
    Run the passed command as a subprocess of the current Python process,
    returning only the exit status of this subprocess.

    This function does _not_ raise exceptions on subprocess failure. To do so,
    consider calling `run()` instead.

    Parameters
    ----------
    command_words : list
        List of one or more shell words comprising this command.
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Returns
    ----------
    int
        Exit status returned by this subprocess.
    '''

    # Avoid circular import dependencies.
    from betse.util.command.exits import FAILURE_DEFAULT

    # Sanitize these arguments.
    _init_run_args(popen_kwargs)

    # Run this command *WITHOUT* raising an exception on command failure.
    try:
        exit_status = subprocess.call(command_words, **popen_kwargs)
    # If this command failed to halt before triggering a timeout, the "timeout"
    # keyword argument was passed *AND* this command has effectively failed.
    # Since the prior call has already guaranteeably terminated this command,
    # this exception is safely convertable into failure exit status.
    except TimeoutExpired:
        exit_status = FAILURE_DEFAULT

    # Return this exit status.
    return exit_status

# ....................{ RUNNERS ~ stdout                   }....................
def run_with_stdout_captured(command_words: list, **popen_kwargs) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    capturing and returning all stdout output by this subprocess _and_ raising
    an exception on subprocess failure.

    Parameters
    ----------
    command_words : list
        List of one or more shell words comprising this command.
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Returns
    ----------
    str
        All stdout captured from this subprocess, stripped of all trailing
        newlines (as under most POSIX shells) _and_ decoded with the current
        locale's preferred encoding (e.g., UTF-8).

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Sanitize these arguments.
    _init_run_args(command_words, popen_kwargs)

    # Capture this command's stdout, raising an exception on command failure
    # (including failure due to an expired timeout).
    command_stdout = subprocess.check_output(command_words, **popen_kwargs)

    # Return this stdout, stripped of all trailing newlines.
    return command_stdout.rstrip('\n')


def run_with_output_interleaved(command_words: list, **popen_kwargs) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    capturing and returning all stdout and stderr output by this subprocess
    (interleaved together) _and_ raising an exception on subprocess failure.

    Parameters
    ----------
    command_words : list
        List of one or more shell words comprising this command.
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Returns
    ----------
    str
        All stdout and stderr captured from this subprocess, interleaved
        together in output order, stripped of all trailing newlines (as under
        most POSIX shells) _and_ decoded with the current locale's preferred
        encoding (e.g., UTF-8).

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Redirect stderr to stdout.
    popen_kwargs['stderr'] = subprocess.STDOUT

    # Capture and return this command's stdout and stderr.
    return run_with_stdout_captured(command_words, **popen_kwargs)

# ....................{ PRIVATE                            }....................
def _init_run_args(command_words: list, popen_kwargs: dict) -> None:
    '''
    Sanitize the dictionary of keyword arguments to be passed to the
    `subprocess.Popen()` callable with sane defaults.

    `close_fds`
    ----------
    If the current platform is vanilla Windows _and_ none of the `stdin`,
    `stdout`, `stderr`, or `close_fds` arguments are passed, the latter argument
    will be explicitly set to `False` -- causing the command to be run to
    inherit all file handles (including stdin, stdout, and stderr) from the
    current process. By default, `subprocess.Popen` documentation insists that:

    > On Windows, if `close_fds` is `True` then no handles will be inherited by
    > the child process.

    The child process will then open new file handles for stdin, stdout, and
    stderr. If the current terminal is a Windows Console, the underlying
    terminal devices and hence file handles will remain the same, in which case
    this is _not_ an issue. If the current terminal is Cygwin-based (e.g.,,
    MinTTY), however, the underlying terminal devices and hence file handles
    will differ, in which case this behaviour prevents interaction between the
    current shell and the vanilla Windows command to be run below. In
    particular, all output from this command will be squelched.

    If at least one of stdin, stdout, or stderr are redirected to a blocking
    pipe, setting `close_fds` to `False` can induce deadlocks under certain
    edge-case scenarios. Since all such file handles default to `None` and hence
    are _not_ redirected in this case, `close_fds` may be safely set to `False`.

    On all other platforms, if `close_fds` is `True`, no file handles _except_
    stdin, stdout, and stderr will be inherited by the child process. This
    function fundamentally differs in subtle (and only slightly documented ways)
    between vanilla Windows and all other platforms. These discrepancies appear
    to be harmful but probably unavoidable, given the philosophical gulf between
    vanilla Windows and all other platforms.

    Arguments
    ----------
    command_words : list
        List of one or more shell words comprising this command.
    popen_kwargs : dict
        Dictionary of keyword arguments to be sanitized.
    '''
    assert types.is_sequence_nonstr_nonempty(command_words), (
        types.assert_not_sequence_nonstr_nonempty(
            command_words, 'Command words'))
    assert types.is_mapping(popen_kwargs), (
        types.assert_not_mapping(popen_kwargs))

    # Avoid circular import dependencies.
    from betse.util.os import oses, shells
    from betse.util.type import dicts

    # Log the command to be run before doing so.
    logs.log_debug('Running command: {}'.format(' '.join(command_words)))

    # If this is vanilla Windows, sanitize the "close_fds" argument.
    if oses.is_windows_vanilla() and not dicts.is_keys(
        popen_kwargs, 'stdin', 'stdout', 'stderr', 'close_fds'):
        popen_kwargs['close_fds'] = False

    # Isolate the current set of environment variables to this command,
    # preventing concurrent changes in these variables in the current process
    # from affecting this command's subprocess.
    popen_kwargs['env'] = shells.get_environment()

    # Decode command output with the current locale's preferred encoding.
    popen_kwargs['universal_newlines'] = True
