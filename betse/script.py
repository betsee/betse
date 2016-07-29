#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# Issue #2: Feature Request - Interactive REPL + Non-interactive Scripting
#
# To facilitate sane non-interactive scripting, the `betse.script` submodule is
# created to make BETSE's various submodules easly accessible to scripts and to
# initialize the BETSE environment.
#
# First perform importations:

# Everyone wants to plot!
import matplotlib as plt
# And everyone wants to NumPy!
import numpy as np          

# For visualization purposes, it is useful to have access to Cells
from betse.science.cells import Cells
# Access to the configuration is very important; Parameters makes this a snap
from betse.science.parameters import Parameters
# And of course scripts will need to access the simulator
from betse.science.sim import Simulator
# It would be great to be able run simulations
from betse.science.simrunner import SimRunner
# Generating logs is useful
import betse.util.io.log.logs as logs
# But typing "logs.<function>" is cluttersome
from betse.util.io.log.logs import log_info, log_exception

# The following includes are generally required for scripts to be useful,
# but as they are undergoing changes, we do not import them (for now).

# This is a must if scripts are going to generate consistent plots
# import betse.science.plot as plot

# Filehandling provides the `loadSim` and `loadInit` functions; necessary
# for... you know... loading simulations and initializeations.
# import betse.science.filehandling as fh

# We also need to initize the BETSE environment. That said, `betse.science`
# already does that at `betse/science/__init__.py:14` by calling
# `ignition.init()`. We include the call below to ensure that it gets called;
# it becomes a noop if it has already been called and makes the intention
# explicit.
from betse import ignition
ignition.init()

# The following functions wrap the standard BETSE API in functions to make
# scripting more convenient and the REPL easier to use.

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
