#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import time
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionSimulation, BetseExceptionParameters
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
    'Running simulation with configuration file "{}".'.format(
        _config_basename))

start_time = time.time()  # get a start value for timing the simulation

p = Parameters(config_filename=_config_filename)  # create an instance of Parameters
p.set_time_profile(p.time_profile_sim)  # force the time profile to be initialize
p.run_sim = True  # set on the fly a boolean to let simulator know we're running a full simulation
sim = Simulator(p)  # create an instance of Simulator

if files.is_file(sim.savedInit):
    sim, cells, p_old = fh.loadSim(sim.savedInit)  # load the initialization from cache
    p.sim_ECM = cells.sim_ECM

    # check to ensure compatibility between original and present sim files:
    if p_old.config['general options'] != p.config['general options'] or\
            p_old.config['world options'] != p.config['world options']:
        # p_old.config['tissue profile definition'] != p.config['tissue profile definition']:
        raise BetseExceptionParameters(
                    'Important config file options are out of sync between '
                    'seed and this init/sim attempt! '
                    'Run "betse seed" again to match the current settings of '
                    'this config file.')

else:
    logs.log_info("No initialization file found to run this simulation!")

    raise BetseExceptionSimulation("Simulation terminated due to missing initialization. Please run "
                                       "an initialization and try again.")

# Reinitialize save and load directories in case params defines new ones
# for this sim.
sim.fileInit(p)

# Run and save the simulation to the cache.
sim.sim_info_report(cells, p)

sim.run_sim_core(cells, p)

logs.log_info(
    'The simulation took {} seconds to complete.'.format(
        round(time.time() - start_time, 2)))
logs.log_info(
    'When ready, close all of the figure windows to end the program.')

if p.turn_all_plots_off is False:
    # As colormaps are deleted from p prior to saving in sim, create a
    # fresh instance of Parameters.
    p = Parameters(config_filename=_config_filename)
    plot_all(cells, sim, p, plot_type='sim')
    plt.show()




