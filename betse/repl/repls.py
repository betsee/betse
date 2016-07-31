#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import betse.pathtree as pathtree
import betse.util.io.log.logs as logs
from betse.exceptions import BetseFunctionException
from betse.repl.environment import repl_env
from betse.util.py import modules

__has_ptpython = modules.is_module('ptpython')

def start_repl(repl_type : str = None):
    '''
    Drop into a REPL.

    This function starts up a REPL. If the `repl_type` argument is `None`, then
    it tries to start `ptpython` and falls back to `code` if this fails. If
    the `repl_type` argument is `not None`, then it attempts to run the
    `repl_type` as specified by the string.

    An exception is raised if `repl_type` is not `None`, `'ptpython'` or`'code'`.

    Parameters
    ----------
    repl_type : str
        The type of REPL to prefer, either `None`, `ptpython` or `code`.
    '''
    if repl_type is None:
        if __has_ptpython:
            start_ptpython_repl()
        else:
            start_code_repl()
    elif repl_type == 'ptpython':
        if __has_ptpython:
            start_ptpython_repl()
        else:
            logs.log_info("The ptpython module does not appear to be installed.")
            logs.log_info("Falling back to a code-based REPL.")
            start_code_repl()
    elif repl_type == 'code':
        start_code_repl()
    else:
        raise BetseFunctionException("unexpected REPL type: \"{}\"".format(repl_type))

def start_code_repl():
    '''
    Start a REPL built around the python `code` module
    '''
    import code
    import readline

    history_filename = pathtree.REPL_HISTORY_FILENAME + ".code"
    readline.set_history_length(1000)
    readline.read_history_file(history_filename)
    try:
        code.interact(banner="", local=repl_env)
    except SystemExit:
        pass
    readline.write_history_file(history_filename)

def start_ptpython_repl():
    '''
    Start a REPL built around the `ptpython` module
    '''
    from ptpython.repl import embed
    try:
        embed(globals=None, locals=repl_env,
            history_filename=pathtree.REPL_HISTORY_FILENAME)
    except SystemExit:
        pass
