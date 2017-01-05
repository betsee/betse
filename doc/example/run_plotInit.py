#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionSimulation
from betse.science import filehandling as fh
from betse.science.parameters import Parameters
from betse.science.plot.pipeline import plot_all
from betse.science.sim import Simulator
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse import ignition, pathtree

# Initialize the current application before doing anything else.
ignition.init()
# at the moment only the default config file can be used
config_filename = pathtree.CONFIG_DEFAULT_FILENAME

# Validate and localize such filename.
files.die_unless_file(config_filename)
_config_filename = config_filename
_config_basename = paths.get_basename(_config_filename)

logs.log_info(
    'Plotting initialization with configuration "{}".'.format(_config_basename))

p = Parameters(config_filename=_config_filename)  # create an instance of Parameters
sim = Simulator(p)  # create an instance of Simulator

if files.is_file(sim.savedInit):
    sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
else:
    raise BetseExceptionSimulation(
        "Ooops! No such initialization file found to plot!")

plot_all(cells, sim, p, plot_type='init')

if p.turn_all_plots_off is False:
    plt.show()