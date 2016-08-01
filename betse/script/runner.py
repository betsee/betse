#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

def run_script(*scripts, dirty = False, globals = globals(), locals = locals()):
    '''
    Run a scripts within the local environment.
    
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
    scripts : *str
        The absolute or relative path of the script, or a sequence of such
        paths.
    
    dirty : bool
        `True` if local changes made in the script should propagate back to
        the caller, or `False` if they should be discarded.

    globals : dict
        A dictionary of variable, function and module bindings to use as the
        global namespace for the script.

    locals : dict
        A dictionary of variables, function and module bindings to use as the
        local namespace for the script. Note that unless `dirty` is `True`,
        any changes made to the local namespace will be discarded.
    '''
    from betse.util.io.log.logs import log_info

    if scripts is ():
        log_info("No scripts to execute")
    else:
        for script in scripts:
            log_info("Executing script: \"{}\"".format(script))
            with open(script) as f:
                if dirty:
                    exec(f.read(), globals, locals)
                else:
                    exec(f.read())
            print()
