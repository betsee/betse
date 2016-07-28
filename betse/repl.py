#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
import code
import betse.util.io.log.logs as logs

def start_repl():
    '''
    Drop into a REPL
    '''
    start_code_repl()

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

    # Create a banner to display as the REPL fires up
    banner = '[betse] Starting the interactive environment'

    # And kick off the REPL. In the event that a `SystemExit` is raised, we
    # catch it and return gracefully. This is because the `quit` function
    # defined above raises such an exception for this purpose exactly.
    try:
        code.interact(banner=banner, readfunc=readfunc, local=env)
    except SystemExit:
        pass
