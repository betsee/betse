#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.
'''
This module provides the environmental context for BETSE REPLs.

Each function and variable in this module is loaded into the `repl_env`
dictionary via a call to `locals`. For safety, `repl_env` is the only
attribute that should be imported from this module.
'''
from betse.script import *

__betse_repl__ = True
'''
Signify that the environment is an active REPL.

It is important that scripts be able to tell that they are running in a REPL
as opposed to standalone mode. For example, if the script's author opts to
follow the "if main" convention, then the script will only run in standalone
mode as `__name__ == "betse.repl.environment"` within the REPL. The
`__betse_repl__` flag makes things a little more concise by allowing:

    if __name__ == '__main__' or __betse_repl__:
        ...

or if behavior should differ between hosted and standalone modes:

    if __name__ == '__main__':
        ...
    elif __betse_repl__:
        ...
'''

def quit():
    '''
    Gracefully exit the REPL, returning control to the caller.
    '''
    raise SystemExit(0)

def run_script(script, *args, dirty = False, g = globals(), l = locals()):
    '''
    Run a script within the local environment, passing *args* as arguments

    If more than one script is provided, each is executed in turn. If any
    script raises an `Exception` then the entire pipeline halts.

    .. caution::
        Note that there is no way to ensure that a script will be perfectly
        "clean" as there can be side-effects beyond those made within the local
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
