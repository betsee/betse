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

def seed(source):
    '''
    Seed a simulation returning a `SimRunner` instance.

    If *source* is an instance of `SimRunner`, then it is used to seed the
    world. Otherwise *source* is first used to initialize a `SimRunner` which
    is then used to seed.

    Parameters
    ----------
    source
        Either an instance of `SimRunner` or something that can initialize
        a `SimRunner`

    Returns
    -------
    An instance of `SimRunner`
    '''
    from betse.util.type import types

    if types.is_simrunner(source):
        source.makeWorld()
        return source
    else:
        runner = SimRunner(source)
        return seed(runner)

def initialize(source):
    '''
    Run an initialization simulation returning a `SimRunner` instance.

    If *source* is an instance of `SimRunner`, then it is used run the
    initialization. Otherwise *source* is first used to initialize a `SimRunner`
    which is then used to run the simulation.

    Parameters
    ----------
    source
        Either an instance of `SimRunner` or something that can initialize
        a `SimRunner`

    Returns
    -------
    An instance of `SimRunner`
    '''
    from betse.util.type import types

    if types.is_simrunner(source):
        source.initialize()
        return source
    else:
        runner = SimRunner(source)
        return init(runner)

def simulate(source):
    '''
    Run a simulation returning a `SimRunner` instance.

    If *source* is an instance of `SimRunner`, then it is used run the
    simulation. Otherwise *source* is first used to initialize a `SimRunner`
    which is then used to run the simulation.

    Parameters
    ----------
    source
        Either an instance of `SimRunner` or something that can initialize
        a `SimRunner`

    Returns
    -------
    An instance of `SimRunner`
    '''
    from betse.util.type import types

    if types.is_simrunner(source):
        source.simulate()
        return source
    else:
        runner = SimRunner(source)
        return sim(runner)

def read_config(config_filename : str):
    '''
    Read a configuration file into a `Parameters` object.

    Parameters
    ----------
    config_filename : str
        The filename of the YAML configuration file

    Returns
    -------
    A `Parameters` instance
    '''
    return Parameters(config_filename)

def load_world(source):
    '''
    Load a world from some *source*.

    Parameters
    ----------
    source
        An instance of `Cells` or `Parameters`, or a path to a YAML
        configuration file.

    Returns
    -------
    (Cells, Parameters)
        A 2-tuple `(sim, cells, p)` as loaded the *source*.
    '''
    from betse.util.type import types
    from betse.science.filehandling import loadWorld

    if types.is_cells(source):
        return loadWorld(source.savedWorld)
    elif types.is_parameters(source):
        return load_world(Cells(source))
    else:
        return load_world(Parameters(source))

def load_init(source):
    '''
    Load an initialization simulation from some *source*.

    Parameters
    ----------
    source
        An instance of `Simulator` or `Parameters`, or a path to a YAML
        configuration file.

    Returns
    -------
    (Simulator, Cells, Parameters)
        A 3-tuple `(sim, cells, p)` as loaded the *source*.
    '''
    from betse.util.type import types
    from betse.science.filehandling import loadSim

    if types.is_simulator(source):
        return loadSim(source.savedInit)
    elif types.is_parameters(source):
        return load_init(Simulator(source))
    else:
        return load_init(Parameters(source))

def load_sim(source):
    '''
    Load a simulation from some *source*.

    Parameters
    ----------
    source
        An instance of `Simulator` or `Parameters`, or a path to a YAML
        configuration file.

    Returns
    -------
    (Simulator, Cells, Parameters)
        A 3-tuple `(sim, cells, p)` as loaded from the *source*.
    '''
    from betse.util.type import types
    from betse.science.filehandling import loadSim

    if types.is_simulator(source):
        return loadSim(source.savedSim)
    elif types.is_parameters(source):
        return load_sim(Simulator(source))
    else:
        return load_sim(Parameters(source))

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
