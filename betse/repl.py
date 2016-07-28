#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import code
import readline
import betse.util.io.log.logs as logs
from betse.exceptions import BetseExceptionFunction

# Attempt to import the ptpython module
try:
    from ptpython.repl import embed
    __has_ptpython = True
except ImportError:
    __has_ptpython = False

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
        raise BetseExceptionFunction("unexpected REPL type: \"{}\"".format(repl_type))

def start_code_repl():
    '''
    Start a REPL built around the python `code` module
    '''
    # We need a function to quit the REPL
    def quit():
        raise SystemExit

    # Create an environment to use as our REPL's namespace, and override the
    # '__name__' variable. This will let use distinguish between interactive
    # and non-interactive environments.
    env = {'__name__' : '__betse_repl__'}

    # Add out local variable and function definitions to the REPL's environment
    env.update(locals())

    # And kick off the REPL. In the event that a `SystemExit` is raised, we
    # catch it and return gracefully. This is because the `quit` function
    # defined above raises such an exception for this purpose exactly.
    try:
        code.interact(banner="", local=env)
    except SystemExit:
        pass

def start_ptpython_repl():
    '''
    Start a REPL built around the `ptpython` module
    '''
    # Create an environment to use as our REPL's namespace, and override the
    # '__name__' variable. This will let use distinguish between interactive
    # and non-interactive environments.
    env = {'__name__' : '__betse_repl__'}

    # And kick off the REPL. In the event that a `SystemExit` is raised, we
    # catch it and return gracefully.
    try:
        embed(globals=globals(), locals=locals())
    except SystemExit:
        pass
