#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level convenience scripting API for end-user purposes.

To facilitate sane non-interactive scripting, this submodule renders BETSE's
various submodules easily accessible to external third-party scripts. For
convenience, BETSE itself is initialized on importing this submodule.
'''

# ....................{ IMPORTS ~ third-party              }....................
# Everyone wants to plot!
import matplotlib as plt
# And everyone wants to NumPy!
import numpy as np

# ....................{ IMPORTS ~ first-party              }....................
# For visualization purposes, it is useful to have access to Cells
from betse.science.cells import Cells
# Access to the configuration is very important; Parameters makes this a snap
from betse.science.parameters import Parameters
# And of course scripts will need to access the simulator
from betse.science.sim import Simulator
# It would be great to be able run simulations
from betse.science.simrunner import SimRunner
# Generating logs is useful
from betse.util.io.log import logs
# But typing "logs.<function>" is cluttersome
from betse.util.io.log.logs import log_info, log_exception

# The following includes are generally required for scripts to be useful,
# but as they are undergoing changes, we do not import them (for now).

# Filehandling provides the `loadSim` and `loadInit` functions; necessary
# for... you know... loading simulations and initializations.
# import betse.science.filehandling as fh

# ....................{ IMPORTS ~ first-party : script     }....................
# The following imports facilitate argument parsing by scripts.
from betse.script.argparse import ArgumentParser, argv

# The following API imports make scripting cleaner and the REPL easier to use.
from betse.script.api import seed, initialize, simulate
# from betse.script.api import read_config, load_world, load_init, load_sim

# ....................{ INITIALIZATIONS                    }....................
# We also need to initialze the BETSE environment. That said, `betse.science`
# already does that at `betse/science/__init__.py:14` by calling
# `ignition.init()`. We include the call below to ensure that it gets called;
# it becomes a noop if it has already been called and makes the intention
# explicit.
from betse import ignition
ignition.init()
