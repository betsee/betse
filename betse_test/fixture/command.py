#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level external command fixtures.

These fixtures automate testing of externally runnable commands in subprocesses
of the current `py.test` process.
'''

# ....................{ IMPORTS                            }....................
from betse.util.py import pys
from betse.util.type import types

# ....................{ RUNNERS                            }....................
#FIXME: Shift all of the following into the main codebase, if not already.

def run_python(command_args: list, **popen_kwargs) -> int:
    '''
    Rerun the currently running Python interpreter as a subprocess of the
    current Python process, returning this subprocess' exit status but neither
    stdout nor stderr.

    Parameters
    ----------
    command_args : list
        List of zero or more arguments to be passed to this interpreter.
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to the `subprocess.Popen()`
        callable internally called by this function. See the `run_command()`
        docstring for further details.

    Returns
    ----------
    int
        Exit status returned by this subprocess.

    See Also
    ----------
    run_command()
        Low-level commentary on subprocess execution.
    '''
    assert types.is_sequence(command_args), (
        types.assert_not_sequence(command_args))

    # List of one or more shell words comprising this command.
    command_line = pys.get_command_line_prefix() + command_args

    #FIXME: Fill me in, please!

    # Rerun this interpreter, propagating its exit status as our own.
    return run_command(command_line, **popen_kwargs)


#FIXME: Implement me, please!

def run_command(command_line: list, **popen_kwargs) -> int:
    '''
    Run the passed command as a subprocess of the current Python process,
    returning this subprocess' exit status but neither stdout nor stderr.

    Parameters
    ----------
    command_line : list
        List of one or more shell words comprising this command, whose:
        * Mandatory first element is either:
            * This command's absolute or relative path.
            * This command's basename, in which case the first command with that
              basename in the current `${PATH}` environment variable will be
              run. If no such command is found, an exception is raised.
        * Optional subsequent elements are this command's arguments (in order).
    popen_kwargs : dict
        Dictionary of keyword arguments to be passed to the `subprocess.Popen()`
        callable internally called by this function. Frequently passed
        arguments of interest include:
        * `cwd`, the absolute path of the current working directory (CWD) from
          which this command is to be run. Defaults to the current CWD.
        * `timeout`, the maximum number of milliseconds this command is to be
          run for. Commands with execution time exceeding this timeout will be
          mercifully killed. Defaults to `None`, in which case this command may
          run indefinitely.

    Returns
    ----------
    int
        Exit status returned by this subprocess.
    '''

    pass
