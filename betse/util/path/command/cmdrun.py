#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Command runners** (i.e., functions running commands on behalf of callers).

Command Words Arguments
----------
Most runners accept a mandatory ``command_words`` parameter, a list of one or
more shell words comprising this command whose:

* Mandatory first item is either:

  * This command's absolute or relative path.
  * This command's basename, in which case the first command with that basename
    in the current ``${PATH}`` environment variable will be run. If no such
    command is found, an exception is raised.

* Optional subsequent items are this command's arguments (in order). Note that
  these arguments are passed as is to a low-level system call rather than
  intprereted by a high-level shell (e.g., ``/bin/sh`` on POSIX-compatible
  platforms) and hence should *not* be shell-quoted. Indeed, shell quoting these
  arguments is likely to result in erroneous command behaviour. The principal
  exception to this heuristic are **GNU-style long value options** (i.e.,
  ``--``-prefixed options accepting a ``=``-delimited value), whose values
  should either be:

  * Passed as a separate argument *without* being shell-quoted.
  * Concatenated to the current ``--``-prefixed argument delimited by ``=`` and
    shell-quoted.

``Popen()`` Keyword Arguments
----------
Most runners accept the same optional keyword arguments accepted by the
:meth:`subprocess.Popen.__init__` constructor, including:

* ``cwd``, the absolute path of the current working directory (CWD) from which
  this command is to be run. Defaults to the current CWD. **Unfortunately, note
  that this keyword argument appears to be erroneously ignored on numerous
  platforms (e.g., Windows XP).** For safety, the
  :func:`betse.util.os.shell import shelldir.setting_cwd` context manager should
  typically be leveraged instead.
* ``timeout``, the maximum number of milliseconds this command is to be run for.
  Commands with execution time exceeding this timeout will be mercifully killed.
  Defaults to ``None``, in which case this command is run indefinitely.
'''

# ....................{ IMPORTS                            }....................
import subprocess
from betse.exceptions import BetseCommandException
from betse.util.io.log import logs
from betse.util.io.log.logenum import LogLevel
from betse.util.type.mapping import mappings
from betse.util.type.types import (
    type_check,
    MappingType,
    MappingOrNoneTypes,
    SequenceTypes,
)
from io import TextIOWrapper
from subprocess import CalledProcessError, Popen, PIPE, TimeoutExpired
from threading import Thread

# ....................{ GLOBALS                            }....................
BUFFER_SIZE_DEFAULT = -1
'''
Default subprocess buffer size for the current platform (synonymous with the
current :data:`io.DEFAULT_BUFFER_SIZE`) suitable for passing as the ``bufsize``
parameter accepted by :meth:`subprocess.Popen.__init__` method.
'''


BUFFER_SIZE_NONE = 0
'''
Unbuffered subprocess buffer size suitable for passing as the ``bufsize``
parameter accepted by :meth:`subprocess.Popen.__init__` method.

Both reading from and writing to an unbuffered subprocess is guaranteed to
perform exactly one system call (``read()`` and ``write()``, respectively) and
can return short (i.e., produce less bytes than the requested number of bytes).
'''


BUFFER_SIZE_LINE = 1
'''
Line-buffered subprocess buffer size suitable for passing as the ``bufsize``
parameter accepted by :meth:`subprocess.Popen.__init__` method.

Reading from a line-buffered subprocess is guaranteed to block until the
subprocess emits a newline, at which point all output emitted between that
newline inclusive and the prior newline exclusive is consumed.
'''

# ....................{ RUNNERS                            }....................
def run_or_die(
    command_words: SequenceTypes, popen_kwargs: MappingOrNoneTypes = None
) -> None:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure *and* forwarding all
    standard output and error output by this subprocess to the standard output
    and error file handles of the current Python process.

    This exception contains the exit status of this subprocess.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import cmdexit

    # Run this command, raising an exception on command failure. For
    # reusability, reimplement the subprocess.check_call() function here rather
    # than explicitly call this function. The latter approach would require
    # duplicating logic between this and the get_exit_status() function.
    exit_status = get_exit_status(
        command_words=command_words, popen_kwargs=popen_kwargs)
    if cmdexit.is_failure(exit_status):
        raise CalledProcessError(exit_status, command_words)

# ....................{ GETTERS                            }....................
def get_exit_status(
    command_words: SequenceTypes, popen_kwargs: MappingOrNoneTypes = None
) -> int:
    '''
    Run the passed command as a subprocess of the current Python process,
    returning only the exit status of this subprocess *and* forwarding all
    standard output and error output by this subprocess to the standard output
    and error file handles of the current Python process.

    This function raises *no* exceptions on subprocess failure. To do so,
    consider calling the :func:`run` function instead.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    Returns
    ----------
    int
        Exit status returned by this subprocess.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command.cmdexit import FAILURE_DEFAULT

    # Sanitize these arguments.
    popen_kwargs = _init_popen_kwargs(command_words, popen_kwargs)

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

# ....................{ GETTERS ~ stdout                   }....................
def get_stdout_or_die(
    command_words: SequenceTypes, popen_kwargs: MappingOrNoneTypes = None
) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure *and* capturing and returning all
    stdout output by this subprocess.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    Returns
    ----------
    str
        All stdout captured from this subprocess, stripped of all trailing
        newlines (as under most POSIX shells) *and* decoded with the current
        locale's preferred encoding (e.g., UTF-8).

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Sanitize these arguments.
    popen_kwargs = _init_popen_kwargs(command_words, popen_kwargs)

    # Capture this command's stdout, raising an exception on command failure
    # (including failure due to an expired timeout).
    command_stdout = subprocess.check_output(command_words, **popen_kwargs)

    # Return this stdout, stripped of all trailing newlines.
    return command_stdout.rstrip('\n')


def get_output_interleaved_or_die(
    command_words: SequenceTypes, popen_kwargs: MappingOrNoneTypes = None
) -> str:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure *and* capturing and returning all
    stdout and stderr output by this subprocess interleaved together (in
    arbitrary order).

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    Returns
    ----------
    str
        All stdout and stderr captured from this subprocess, interleaved
        together in output order, stripped of all trailing newlines (as under
        most POSIX shells) *and* decoded with the current locale's preferred
        encoding (e.g., UTF-8).

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # If these keyword arguments are empty, default to the empty dictionary.
    if not popen_kwargs:
        popen_kwargs = {}

    # Redirect stderr to stdout.
    popen_kwargs['stderr'] = subprocess.STDOUT

    # Capture and return this command's stdout and stderr.
    return get_stdout_or_die(
        command_words=command_words, popen_kwargs=popen_kwargs)

# ....................{ LOGGERS                            }....................
def log_output_or_die(
    # Mandatory arguments.
    command_words: SequenceTypes,

    # Optional arguments.
    stdout_log_level: LogLevel = LogLevel.INFO,
    stderr_log_level: LogLevel = LogLevel.ERROR,
    popen_kwargs: MappingOrNoneTypes = None,
) -> None:
    '''
    Run the passed command as a subprocess of the current Python process,
    raising an exception on subprocess failure *and* capturing and logging all
    stdout output by this subprocess with logging level ``INFO`` and all
    stderr output by this subprocess with logging level ``ERROR``.

    For both space efficiency and logging responsiveness, this function captures
    and logs subprocess output in a line-buffered manner. Specifically:

    * Each line of stdout output by this subprocess is logged with logging level
      :attr:`LogLevel.INFO` with multi-threaded asynchronous I/O.
    * Each line of stdout output by this subprocess is logged with logging level
      :attr:`LogLevel.INFO` with multi-threaded asynchronous I/O.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    stdout_log_level : LogLevel
        Logging level with which all stdout output is logged. Defaults to
        :attr:`LogLevel.INFO`.
    stderr_log_level : LogLevel
        Logging level with which all stderr output is logged. Defaults to
        :attr:`LogLevel.ERROR`.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method. Defaults to ``None``, in which
        case the empty dictionary is assumed.

    Raises
    ----------
    CalledProcessError
        Exception raised on subprocess failure.
    '''

    # Sanitize these arguments.
    popen_kwargs = _init_popen_kwargs(command_words, popen_kwargs)

    # Subprocess forked from this process, redirecting both stdout and stderr to
    # consumable pipes in a line-buffered manner.
    command_subprocess = Popen(
        args=command_words,
        stdout=PIPE,
        stderr=PIPE,
        bufsize=BUFFER_SIZE_LINE,
        **popen_kwargs
    )

    # Threads logging each line of stdout and stderr output by this subprocess
    # with the appropriate logging levels.
    stdout_logger = Thread(
        target=_log_pipe_lines,
        args=(command_subprocess.stdout, stdout_log_level))
    stderr_logger = Thread(
        target=_log_pipe_lines,
        args=(command_subprocess.stderr, stderr_log_level))

    # Run this subprocess in a multi-threaded manner logged by these
    # asynchronous threads, raising an exception on subprocess failure.
    stdout_logger.start()
    stderr_logger.start()


@type_check
def _log_pipe_lines(pipe: TextIOWrapper, log_level: LogLevel) -> None:
    '''
    Iteratively log each line buffered by the passed subprocess pipe as a new
    log message with the passed log level.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : Mapping
        Dictionary of keyword arguments to pass to the
        :meth:`subprocess.Popen.__init__` method.
    '''

    # Avoid circular import dependencies.
    from betse.util.io.iofiles import READLINE_EOF
    from betse.util.type.text import strs

    # With this pipe contextually opened for reading...
    with pipe:
        # For each plaintext line emitted by this pipe (falling back to the
        # standard EOF emitted by the readline() method in the event of
        # unexpected pipe closure)...
        for line in iter(pipe.readline, READLINE_EOF):
            # If this line signifies end-of-file (EOF), stop consuming lines.
            if line == READLINE_EOF:
                break

            # Message to be logged, stripped of trailing newline if any.
            message = strs.remove_newlines_suffix(line)

            # Log this line as a new log message with this logging level.
            logs.log_levelled(message, log_level)

# ....................{ PRIVATE                            }....................
@type_check
def _init_popen_kwargs(
    command_words: SequenceTypes, popen_kwargs: MappingOrNoneTypes
) -> MappingType:
    '''
    Sanitized dictionary of all keyword arguments to pass to the
    :class:`subprocess.Popen` callable when running the command specified by the
    passed shell words with the passed user-defined keyword arguments.

    `close_fds`
    ----------
    If the current platform is vanilla Windows *and* none of the ``stdin``,
    ``stdout``, ``stderr``, or ``close_fds`` arguments are passed, the latter
    argument will be explicitly set to ``False`` -- causing the command to be
    run to inherit all file handles (including stdin, stdout, and stderr) from
    the current process. By default, :class:`subprocess.Popen` documentation
    insists that:

    > On Windows, if ``close_fds`` is ``True`` then no handles will be inherited
    > by the child process.

    The child process will then open new file handles for stdin, stdout, and
    stderr. If the current terminal is a Windows Console, the underlying
    terminal devices and hence file handles will remain the same, in which case
    this is *not* an issue. If the current terminal is Cygwin-based (e.g.,,
    MinTTY), however, the underlying terminal devices and hence file handles
    will differ, in which case this behaviour prevents interaction between the
    current shell and the vanilla Windows command to be run below. In
    particular, all output from this command will be squelched.

    If at least one of stdin, stdout, or stderr are redirected to a blocking
    pipe, setting ``close_fds`` to ``False`` can induce deadlocks under certain
    edge-case scenarios. Since all such file handles default to ``None`` and
    hence are *not* redirected in this case, ``close_fds`` may be safely set to
    ``False``.

    On all other platforms, if ``close_fds`` is ``True``, no file handles
    *except* stdin, stdout, and stderr will be inherited by the child process.
    This function fundamentally differs in subtle (and only slightly documented
    ways) between vanilla Windows and all other platforms. These discrepancies
    appear to be harmful but probably unavoidable, given the philosophical gulf
    between vanilla Windows and all other platforms.

    Parameters
    ----------
    command_words : SequenceTypes
        List of one or more shell words comprising this command.
    popen_kwargs : optional[MappingType]
        Dictionary of all keyword arguments to be sanitized if any *or* ``None``
        otherwise, in which case the empty dictionary is defaulted to.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import cmds
    from betse.util.os import oses
    from betse.util.os.shell import shellenv

    # If this list of shell words is empty, raise an exception.
    if not command_words:
        raise BetseCommandException('Non-empty command expected.')

    # If these keyword arguments are empty, default to the empty dictionary.
    if popen_kwargs is None:
        popen_kwargs = {}

    # If the first shell word is this list is unrunnable, raise an exception.
    cmds.die_unless_command(command_words[0])

    # Log the command to be run before doing so.
    logs.log_debug('Running command: %s', ' '.join(command_words))

    # If this is vanilla Windows, sanitize the "close_fds" argument.
    if oses.is_windows_vanilla() and not mappings.is_keys(
        popen_kwargs, 'stdin', 'stdout', 'stderr', 'close_fds'):
        popen_kwargs['close_fds'] = False

    # Isolate the current set of environment variables to this command,
    # preventing concurrent changes in these variables in the current process
    # from affecting this command's subprocess.
    popen_kwargs['env'] = shellenv.get_env()

    # Decode command output with the current locale's preferred encoding.
    popen_kwargs['universal_newlines'] = True

    # Return these keyword arguments.
    return popen_kwargs
