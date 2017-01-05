#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
**Command runners** (i.e., functions running commands on behalf of callers).

Command Words Arguments
----------
Most runners accept a mandatory `command_words` parameter, a list of one or more
shell words comprising this command whose:

* Mandatory first element is either:
  * This command's absolute or relative path.
  * This command's basename, in which case the first command with that basename
    in the current `${PATH}` environment variable will be run. If no such
    command is found, an exception is raised.
* Optional subsequent elements are this command's arguments (in order).

`Popen()` Keyword Arguments
----------
Most runners accept the same optional keyword arguments accepted by the
`subprocess.Popen.__init__()` constructor, including:

* `cwd`, the absolute path of the current working directory (CWD) from which
  this command is to be run. Defaults to the current CWD.
* `timeout`, the maximum number of milliseconds this command is to be run for.
  Commands with execution time exceeding this timeout will be mercifully killed.
  Defaults to `None`, in which case this command may run indefinitely.
'''

# ....................{ IMPORTS                            }....................
import subprocess
from betse.exceptions import BetseCommandException
from betse.util.io.log import logs
from betse.util.type.types import type_check, MappingType, SequenceTypes
from subprocess import CalledProcessError, TimeoutExpired

# ....................{ RUNNERS                            }....................
def run(command_words: SequenceTypes, **popen_kwargs) -> None:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure.

    This exception contains the exit status of this subprocess.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : Mapping
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import exits

    # Run this command, raising an exception on command failure. For
    # reusability, reimplement the subprocess.check_call() function here rather
    # than explicitly call this function. The latter approach would require
    # duplicating logic between this and the run_nonfatal() function.
    exit_status = run_nonfatal(command_words, popen_kwargs)
    if exits.is_failure(exit_status):
        raise CalledProcessError(exit_status, command_words)


def run_nonfatal(command_words: SequenceTypes, **popen_kwargs) -> int:
    '''
    Run the passed command as a subprocess of the current Python process,
    returning only the exit status of this subprocess.

    This function does _not_ raise exceptions on subprocess failure. To do so,
    consider calling `run()` instead.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : Mapping
        Dictionary of keyword arguments to be passed to `subprocess.Popen()`.

    Returns
    ----------
    int
        Exit status returned by this subprocess.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import FAILURE_DEFAULT

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
def run_capturing_stdout(
    command_words: SequenceTypes, **popen_kwargs) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    capturing and returning all stdout output by this subprocess _and_ raising
    an exception on subprocess failure.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : Mapping
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


def run_interleaving_output(
    command_words: SequenceTypes, **popen_kwargs) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    capturing and returning all stdout and stderr output by this subprocess
    interleaved together (in arbitrary order) _and_ raising an exception on
    subprocess failure.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : Mapping
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
    return run_capturing_stdout(command_words, **popen_kwargs)

# ....................{ PRIVATE                            }....................
@type_check
def _init_run_args(
    command_words: SequenceTypes, popen_kwargs: MappingType) -> None:
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

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : MappingType
        Dictionary of keyword arguments to be sanitized.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import commands
    from betse.util.os import oses
    from betse.util.os.shell import envs
    from betse.util.type import mappings

    # If the passed list of shell words is empty, raise an exception.
    if not command_words:
        raise BetseCommandException('Non-empty command expected.')

    # If the first shell word is this list is unrunnable, raise an exception.
    commands.die_unless_command(command_words[0])

    # Log the command to be run before doing so.
    logs.log_debug('Running command: %s', ' '.join(command_words))

    # If this is vanilla Windows, sanitize the "close_fds" argument.
    if oses.is_windows_vanilla() and not mappings.is_keys(
        popen_kwargs, 'stdin', 'stdout', 'stderr', 'close_fds'):
        popen_kwargs['close_fds'] = False

    # Isolate the current set of environment variables to this command,
    # preventing concurrent changes in these variables in the current process
    # from affecting this command's subprocess.
    popen_kwargs['env'] = envs.get_env()

    # Decode command output with the current locale's preferred encoding.
    popen_kwargs['universal_newlines'] = True
