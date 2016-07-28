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
import matplotlib.pyplot as plt
import numpy as np

from betse.science.simrunner import SimRunner
from betse.science.parameters import Parameters
from betse.science.sim import Simulator

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
    if isinstance(source, SimRunner):
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
    if isinstance(source, SimRunner):
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
    if isinstance(source, SimRunner):
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

repl_env = locals()
