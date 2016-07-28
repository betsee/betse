#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import code
import betse.util.io.log.logs as logs
import betse.exceptions as ex

# Attempt to import the ptpython module
try:
    import ptpython
    __repl_type = 'ptpython'
except ImportError:
    __repl_type = 'code'

def start_repl():
    '''
    Drop into a REPL
    '''
    if __repl_type is 'ptpython':
        start_ptpython_repl()
    elif __repl_type is 'code':
        start_code_repl()
    else:
        raise ex.BetseExceptionModule("\"__repl_type\" has an unexpected value: {}".format(__repl_type))

def start_code_repl():
    '''
    Start a REPL built around the python `code` module
    '''
    # We need a function to quit the REPL
    def quit():
        raise SystemExit

    # And a function to provide a nice BETSE-specific prompt
    def readfunc(prompt):
        return input('[betse-repl] > ')

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
        code.interact(banner="", readfunc=readfunc, local=env)
    except SystemExit:
        pass

def start_ptpython_repl():
    '''
    Start a REPL built around the `ptpython` module
    '''
    logs.log_info("The ptpython-based repl is not yet implemented.")
    logs.log_info("Falling back to a python code-based repl.")
    start_code_repl()
