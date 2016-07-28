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

def init(source):
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

def sim(source):
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

repl_env = locals()
