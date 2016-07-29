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
    raise SystemExit

def run_script(scripts, dirty=False, globals=globals(), locals=locals()):
    '''
    Run a script (or sequence of scripts) within the local environment.
    
    If `scripts` is a sequence of scripts, each is executed in turn. If any
    script raises an `Exception` then the entire pipeline halts.

    .. caution::
        Note that there is no way to ensure that a script will be perfectly
        "clean" as there can be side-effects beyond those make within the local
        python environment. For example, the script could write to one of the
        various simulation files. The author sees no good way to ensure that
        this does not happen.

    Parameters
    ----------
    scripts : str or sequence
        The absolute or relative path of the script, or a sequence of such
        paths.
    
    dirty : bool
        `True` if local changes made in the script should propagate back to
        the caller, or `False` if they should be discarded.

    globals : dict
        A dictionary of variable, function and module bindings to use as the
        global namespace for the script.

    locals : dict
        A dictionary of variable, function and module bindings to use as the
        local namespace for the script. Note that unless `dirty` is `True`,
        any changes made to the local namespace will be discarded.
    '''
    from betse.util.type import types
    from betse.exceptions import BetseFunctionException

    def run_single_script(filename):
        with open(filename) as f:
            if dirty:
                exec(f.read(), globals, locals)
            else:
                exec(f.read())

    if types.is_sequence_nonstr(scripts):
        for script in scripts:
            log_info("Executing script: \"{}\"".format(script))
            run_single_script(script)
            print()

    elif types.is_str(scripts):
        run_single_script(scripts)

    else:
        msg = "expected a string or sequence of strings as an argument"
        raise BetseExceptionFunction(msg)

repl_env = locals()
