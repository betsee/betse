#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Dependency-specific REPL facilities.

This module defines functions specific to both official and optional third-party
REPL packages (e.g., :mod:`code`, :mod:`ptpython`).
'''

# ....................{ IMPORTS                            }....................
from betse import pathtree
from betse.cli.repl import environment
from betse.lib import libs
from betse.util.io.log import logs
from betse.util.path import files
from betse.util.path.command import cmdexit
from betse.util.type.types import type_check
from enum import Enum

# ....................{ ENUMS                              }....................
# One-liners are happy liners.
REPLType = Enum('REPLType', ('first_available', 'ptpython', 'code',))
'''
Enumeration of all possible REPLs currently supported by this submodule.
'''

# ....................{ FUNCTIONS                          }....................
@type_check
def start_repl(repl_type: REPLType = REPLType.first_available) -> None:
    '''
    Start a REPL of the passed type.

    Parameters
    ----------
    repl_type : optional[REPLType]
        Type of REPL to prefer. If :data:`REPLType.first_available`, the set of
        all possible REPLs is iteratively searched for the first available REPL;
        else if this REPL is unavailable, the first available REPL is used.
        Defaults to :data:`REPLType.first_available`.
    '''

    if repl_type is REPLType.first_available:
        start_first_repl()
    elif repl_type is REPLType.ptpython:
        start_ptpython_repl()
    elif repl_type is REPLType.code:
        start_code_repl()
    else:
        logs.log_warning(
            'REPL type "{}" unrecognized. '
            'Deferring to first available REPL.'.format(repl_type.name))
        start_first_repl()


def start_first_repl() -> None:
    '''
    Start the first available REPL.
    '''

    if libs.is_runtime_optional('ptpython'):
        start_ptpython_repl()
    else:
        start_code_repl()


def start_ptpython_repl() -> None:
    '''
    Start a REPL based on the optional third-party :mod:`ptpython` package.

    If this package is unavailable, this function defers to the first available
    REPL with a non-fatal warning.
    '''

    # If "ptpython" is unavailable...
    if not libs.is_runtime_optional('ptpython'):
        # Log a non-fatal warning.
        logs.log_warning(
            '"ptpython" module not found. Deferring to first available REPL.')

        # Defer to the first available REPL.
        start_first_repl()

        # Get us out of here, Great Captain.
        return
    # Else, "ptpython" is available.

    # Log this invocation.
    logs.log_info('Starting "ptpython"-based REPL...')

    # Defer heavyweight imports.
    from ptpython.repl import embed

    # If the "ptpython" key is missing from the dictionary of history
    # filenames, then default to no history file. This prevents the readline
    # history files being corrupted by ptpython's unique format.
    history_filename = pathtree.get_repl_history_filename('ptpython')

    # Run this REPL.
    try:
        embed(
            globals=None,
            locals=environment.repl_env,
            history_filename=history_filename,
        )
    # When this REPL halts with error, reraise this exception.
    except SystemExit as exit:
        if cmdexit.is_failure(exit.code):
            raise


def start_code_repl() -> None:
    '''
    Start a REPL based on the canonical :mod:`code` module.
    '''

    # Log this invocation, including the same helpful help line preceding the
    # official Python 3 REPL.
    logs.log_info('Starting "code"-based REPL...')
    logs.log_info(
        'Type "help", "copyright", "credits" or '
        '"license" for more information.')

    # Defer heavyweight imports.
    import code, readline

    # Absolute path of the file persisting a REPL-specific history of commands.
    # Note this REPL leverages a "readline"-style history file format.
    history_filename = pathtree.get_repl_history_filename('readline')
    readline.set_history_length(1000)

    # If this file exists, deserialize this REPL's history from this file.
    if files.is_file(history_filename):
        logs.log_debug('Restoring REPL history from: %s', history_filename)
        readline.read_history_file(history_filename)

    # Run this REPL.
    try:
        code.interact(banner="", local=environment.repl_env)
    # When this REPL halts...
    except SystemExit as exit:
        # If this REPL halted with error, reraise this exception.
        if cmdexit.is_failure(exit.code):
            raise
        # Else, this REPL halted without error. Silently ignore this exception.
    # Serialize this REPL's history back to disk regardless of whether an
    # exception was raised.
    finally:
        logs.log_debug('Preserving REPL history to: %s', history_filename)
        readline.write_history_file(history_filename)
