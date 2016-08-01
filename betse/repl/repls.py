#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import betse.pathtree as pathtree
from betse.repl.environment import repl_env
from betse.util.io.log.logs import log_info, log_warning
from betse.util.py import modules
from betse.util.type.types import type_check
from enum import Enum

__has_ptpython = modules.is_module('ptpython')

class REPLType(Enum):
    '''
    An enumeration of possible REPLs available for use.
    '''
    first_available = 0 # Choose the first available REPL
    ptpython = 1        # Prefer a ptpython REPL
    code = 2            # Prefer a code-based REPL (always available)

@type_check
def start_repl(repl_type : REPLType = REPLType.first_available) -> None:
    '''
    Drop into a REPL.

    This function starts up a REPL. If the `repl_type` is
    `REPLType.first_available`, then each possible REPL is tried in turn
    until one can be successfully started. If the preferred REPL is
    unavailable, then the first available is used.

    Parameters
    ----------
    repl_type : REPLType
        The type of REPL to prefer
    '''
    if repl_type is REPLType.first_available:
        start_first_repl()
    elif repl_type is REPLType.ptpython:
        start_ptpython_repl()
    elif repl_type is REPLType.code:
        start_code_repl()
    else:
        log_warning("The \"{}\" REPL type is not yet implemented.".format(repl_type.name))
        log_warning("Falling back to the first available REPL.")
        start_first_repl()

def start_first_repl():
    '''
    Start the first available REPL.
    '''
    if __has_ptpython:
        start_ptpython_repl()
    else:
        start_code_repl()

def start_ptpython_repl():
    '''
    Start a REPL built around the `ptpython` module.

    If the `ptpython` module is unavailable, then fall back to the first
    available REPL.
    '''
    if not __has_ptpython:
        log_warning("The ptpython module does not appear to be installed.")
        log_warning("Falling back to the first available REPL.")
        start_first_repl()
    else:
        log_info("Starting a ptpython-based REPL.")
        from ptpython.repl import embed
        try:
            embed(globals=None, locals=repl_env,
                history_filename=pathtree.REPL_HISTORY_FILENAME)
        except SystemExit as exit:
            from betse.util.path.command import exits
            if exits.is_failure(exit.code):
                raise

def start_code_repl():
    '''
    Start a REPL built around the python `code` module.
    '''
    log_info("Starting a code-based REPL.")

    import code
    import readline

    history_filename = pathtree.REPL_HISTORY_FILENAME + ".code"
    readline.set_history_length(1000)
    readline.read_history_file(history_filename)
    try:
        code.interact(banner="", local=repl_env)
    except SystemExit as exit:
        from betse.util.path.command import exits
        if exits.is_failure(exit.code):
            raise
    readline.write_history_file(history_filename)
