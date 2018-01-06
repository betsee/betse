#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
from betse.science.cells import Cells
from betse.science.filehandling import loadWorld, loadSim
from betse.science.parameters import Parameters
from betse.science.sim import Simulator
from betse.science.simrunner import SimRunner
from betse.util.type import types
from betse.util.type.types import type_check

# ....................{ PHASES                             }....................
@type_check
def seed(source: (SimRunner, str)) -> SimRunner:
    '''
    Seed a simulation and return the `SimRunner` instance used to do so.

    If `source` is a `SimRunner` instance, this is used to seed; else, `source`
    _must_ be a path to a YAML-formatted simulation configuration file used to
    instantiate a `SimRunner` instance then used to seed.

    Parameters
    ----------
    source : (SimRunner, str)
        Either:
        * An instance of `SimRunner`.
        * The absolute or relative path of a YAML-formatted simulation
          configuration file.

    Returns
    -------
    SimRunner
        Simulation runner running this phase.
    '''

    if types.is_simrunner(source):
        source.make_world()
        return source
    else:
        runner = SimRunner(conf_filename=source)
        return seed(runner)

@type_check
def initialize(source : (SimRunner, str)) -> SimRunner:
    '''
    Run an initialization simulation returning a `SimRunner` instance.

    If *source* is an instance of `SimRunner`, then it is used run the
    initialization. Otherwise *source* is expected to be a path to a
    YAML configuration file, and is used to initialize a `SimRunner`
    which is then used to run the simulation.

    Parameters
    ----------
    source
        Either an instance of `SimRunner` or the path to a YAML configuration
        file

    Returns
    -------
    An instance of `SimRunner`
    '''
    if types.is_simrunner(source):
        source.init()
        return source
    else:
        runner = SimRunner(source)
        return initialize(runner)

@type_check
def simulate(source : (SimRunner, str)) -> SimRunner:
    '''
    Run a simulation returning a `SimRunner` instance.

    If *source* is an instance of `SimRunner`, then it is used run the
    simulation. Otherwise *source* is expected to be a path to a YAML
    configuration file, and is to initialize a `SimRunner` which is then used
    to run the simulation.

    Parameters
    ----------
    source
        Either an instance of `SimRunner` or the path to a YAML configuration
        file

    Returns
    -------
    An instance of `SimRunner`
    '''
    if types.is_simrunner(source):
        source.sim()
        return source
    else:
        runner = SimRunner(source)
        return simulate(runner)

#FIXME: The following functionality no longer complies with the existing BETSE
#API -- particularly, the "Parameters" API. Since the "betse.script" subpackage
#is only currently imported by the "betse.cli.repl" subpackage, since no one
#actually appears to use the "betse repl" enviroonment, *AND* since we currently
#lack sufficient resources to maintain this logic, we should now begin the
#laborous process of excising all of the following subpackages from the
#codebase:
#
#* "betse.script".
#* "betse.cli.repl".

# @type_check
# def read_config(conf_filename : str) -> Parameters:
#     '''
#     Read a configuration file into a `Parameters` object.
#
#     Parameters
#     ----------
#     conf_filename : str
#         The filename of the YAML configuration file
#
#     Returns
#     -------
#     A `Parameters` instance
#     '''
#     return Parameters(conf_filename)
#
# @type_check
# def load_world(source : (Cells, Parameters, str)) -> tuple:
#     '''
#     Load a world from some *source*.
#
#     Parameters
#     ----------
#     source
#         An instance of `Cells` or `Parameters`, or a path to a YAML
#         configuration file
#
#     Returns
#     -------
#     (Cells, Parameters)
#         A 2-tuple `(cells, p)` as loaded from the *source*
#     '''
#     if types.is_cells(source):
#         return loadWorld(source.savedWorld)
#     elif types.is_parameters(source):
#         return load_world(Cells(source))
#     else:
#         return load_world(Parameters(source))
#
# @type_check
# def load_init(source : (Simulator, Parameters, str)) -> tuple:
#     '''
#     Load an initialization simulation from some *source*.
#
#     Parameters
#     ----------
#     source
#         An instance of `Simulator` or `Parameters`, or a path to a YAML
#         configuration file
#
#     Returns
#     -------
#     (Simulator, Cells, Parameters)
#         A 3-tuple `(sim, cells, p)` as loaded from the *source*
#     '''
#     if types.is_simulator(source):
#         return loadSim(source.savedInit)
#     elif types.is_parameters(source):
#         return load_init(Simulator(source))
#     else:
#         return load_init(Parameters(source))
#
# @type_check
# def load_sim(source : (Simulator, Parameters, str)) -> tuple:
#     '''
#     Load a simulation from some *source*.
#
#     Parameters
#     ----------
#     source
#         An instance of `Simulator` or `Parameters`, or a path to a YAML
#         configuration file
#
#     Returns
#     -------
#     (Simulator, Cells, Parameters)
#         A 3-tuple `(sim, cells, p)` as loaded from the *source*
#     '''
#     if types.is_simulator(source):
#         return loadSim(source.savedSim)
#     elif types.is_parameters(source):
#         return load_sim(Simulator(source))
#     else:
#         return load_sim(Parameters(source))
