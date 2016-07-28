#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import code
import readline
import betse.util.io.log.logs as logs

# Attempt to import the ptpython module
try:
    from ptpython.repl import embed
    __has_ptpython = True
except ImportError:
    __has_ptpython = False

def start_repl():
    '''
    Drop into a REPL
    '''
    if __has_ptpython:
        start_ptpython_repl()
    else:
        start_code_repl()

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
