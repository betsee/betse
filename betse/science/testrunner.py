#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

from betse.science.simrunner import SimRunner

# If this module is imported from the command line, run; else, noop.
if __name__ == '__main__':
    # Simulation specified by the default configuration file.
    # sim_runner = SimRunner(
    #     conf_filename = pathtree.get_sim_config_default_filename())

    # print(pathtree.get_sim_config_default_filename())

    sim_runner = SimRunner(
        conf_filename = '/home/pietakio/Documents/BETSE_Study/GRN/test_gene.yaml')

    # Run such simulation.
    # sim_runner.initialize()
    sim_runner.sim_grn()
    # sim_runner.plot_grn()










