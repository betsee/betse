#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

from betse import ignition, pathtree
from betse.science.simrunner import SimRunner

# If this module is imported from the command line, run; else, noop.
if __name__ == '__main__':
    # Initialize the current application *BEFORE* doing anything else.
    # ignition.init()

    # Simulation specified by the default configuration file.
    # sim_runner = SimRunner(
    #     config_filename = pathtree.CONFIG_DEFAULT_FILENAME)

    # print(pathtree.CONFIG_DEFAULT_FILENAME)

    sim_runner = SimRunner(
        config_filename = '/home/pietakio/Documents/BETSE_Study/Metabolism/test_metabo.yaml')

    # Run such simulation.
    # sim_runner.initialize()
    sim_runner.sim_reactions()
    # sim_runner.plot_reactions()










