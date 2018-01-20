#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Environmental context for the BETSE REPL.

Each function and variable in this module is loaded into the global
:data:`repl_env` dictionary via a call to :func:`locals`. For safety, this
dictionary is the only attribute that should be imported from this module.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check

# ....................{ GLOBALS                            }....................
__betse_repl__ = True
'''
`True` only if the current environment is an active REPL.

It is important that scripts be able to tell that they are running in a REPL
as opposed to standalone mode. For example, if the script's author opts to
follow the "if main" convention, then the script will only run in standalone
mode as `__name__ == "betse.repl.environment"` within the REPL. This boolean
simplifies matters by allowing logic resembling:

    if __name__ == '__main__' or __betse_repl__:
        ...

or if behavior should differ between hosted and standalone modes:

    if __name__ == '__main__':
        ...
    elif __betse_repl__:
        ...
'''

# ....................{ FUNCTIONS                          }....................
def quit() -> None:
    '''
    Gracefully exit the current REPL, returning control to the caller.
    '''

    # Avoid polluting the module namespace captured into the "repl_env" global
    # below with these function-specific attributes.
    from betse.util.path.command import cmdexit

    # Exit us up.
    raise cmdexit.raise_success()


#FIXME: Single-letter parameters cause me to poke my eyes out. Would it be
#feasible to rename these to "globals" and "locals" or... well, pretty much
#anything other than "g" and "l"?
@type_check
def run_script(
    script: str,
    *args: str,
    dirty: bool = False,
    g: dict = None,
    l: dict = None
) -> None:
    '''
    Run a script within the local environment passed the passed arguments.

    If more than one script is provided, each is executed in turn. If any
    script raises an exception then the entire pipeline halts.

    .. caution::
        Note that there is no way to ensure that a script will be perfectly
        "clean" as there can be side-effects beyond those made within the local
        Python environment. For example, the script could write to one of the
        various simulation files. Cleanliness is, ultimately, a responsibility
        of the caller.

    Parameters
    ----------
    script : str
        Absolute or relative path of the script to be run.
    args : tuple
        Argument list to be passed to this script.
    dirty : optional[bool]
        If:
        * `True`, local changes made by this script will propagate back to the
          caller _after_ this script halts.
        * `False`, such changes will be discarded.
        Defaults to `False`.
    g : optional[dict]
        Dictionary of variable, function and module bindings to use as the
        global namespace for this script. Defaults to `None`, in which case the
        current :func:`globals` dictionary will be used.
    l : optional[dict]
        Dictionary of variables, function and module bindings to use as the
        local namespace for this script. Note that unless `dirty` is `True`,
        any changes made to this namespace will be discarded. Defaults to
        `None`, in which case the current :data:`repl_env` dictionary will be
        used.
    '''

    # Avoid polluting the module namespace captured into the "repl_env" global
    # below with these function-specific attributes.
    from betse.exceptions import BetseCLIArgParserException
    from betse.script import argv
    from betse.util.io.log.logs import log_info

    # Log this script execution.
    log_info('Running script: %s', script)

    # Record the passed argument list for subsequent retrieval by scripts.
    argv.set_args(script, *args)

    # Run this script.
    try:
        with open(script) as script_file:
            if dirty:
                # Default the passed globals and locals if needed. Due to the
                # default parameter trap, these defaults must *NOT* be defined
                # in this function's signature.
                if g is None:
                    g = globals()
                if l is None:
                    l = repl_env

                # Run this script with access to these globals and locals.
                exec(script_file.read(), g, l)
            else:
                exec(script_file.read())
    except BetseCLIArgParserException:
        pass
    finally:
        argv.uninitialize()
        print()

# ....................{ GLOBALS ~ env                      }....................
# Defer importing all attributes defined by the scripting API until immediately
# *BEFORE* capturing these attributes into the "repl_env" global. Why? Because
# the star-import prevents sane code analysis (e.g., by IDEs) and hence should
# be deferred as long as feasible.
from betse.script import *

repl_env = locals()
'''
Dictionary mapping from the names to values of all attributes defined by this
submodule, including those defined by the :mod:`betse.script` subpackage.
'''
