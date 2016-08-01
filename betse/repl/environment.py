#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
'''
This module provides the environmental context for BETSE REPLs.

Each function and variable in this module is loaded into the `repl_env`
dictionary via a call to `locals`. This is the only symbol that should
be imported from this module.
'''
from betse.script import *

__betse_repl__ = True

def quit():
    '''
    Gracefully exit the REPL, returning control the the caller.
    '''
    raise SystemExit(0)

def run_script(script, *args, dirty = False, g = globals(), l = locals()):
    '''
    Run a script within the local environment, passing *args* as arguments

    If more than one script is provided, each is executed in turn. If any
    script raises an `Exception` then the entire pipeline halts.

    .. caution::
        Note that there is no way to ensure that a script will be perfectly
        "clean" as there can be side-effects beyond those make within the local
        python environment. For example, the script could write to one of the
        various simulation files. The author sees no good way to ensure that
        this does not happen.

    Parameters
    ----------
    script : str
        The absolute or relative path of the script

    args : *str
        The arguments to pass to the script

    dirty : bool
        `True` if local changes made in the script should propagate back to
        the caller, or `False` if they should be discarded.

    g : dict
        A dictionary of variable, function and module bindings to use as the
        global namespace for the script.

    l : dict
        A dictionary of variables, function and module bindings to use as the
        local namespace for the script. Note that unless `dirty` is `True`,
        any changes made to the local namespace will be discarded.
    '''
    from betse.exceptions import BetseArgumentParserException
    from betse.script import argv
    from betse.util.io.log.logs import log_info

    log_info("Executing script \"{}\"".format(script))

    argv.set_args(script, *args)
    try:
        with open(script) as f:
            if dirty:
                exec(f.read(), g, l)
            else:
                exec(f.read())
    except BetseArgumentParserException:
        pass
    finally:
        argv.uninitialize()
        print()

repl_env = locals()
