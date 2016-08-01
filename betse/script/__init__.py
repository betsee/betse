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

# The following imports facilitate argument parsing by scripts
from .argparse import ArgumentParser, betse_argv

# The following API imports make scripting cleaner and the REPL easier to use
from .api import seed, initialize, simulate
from .api import read_config, load_world, load_init, load_sim
