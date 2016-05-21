#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ TESTS                              }....................
#FIXME: While useful, this should be split apart into its constituent steps,
#most of which should probably become fixtures for efficient reuse. Wait, no.
#Ideally, we should construct a linear chain of tests, each exercising a
#simulation phase (e.g., "seed", "init", "sim") depending on the success of the
#prior test in this chain.

def test_cli_sim_default(betse_cli, betse_sim_config_default) -> None:
    '''
    Test the simulation of the default simulation configuration.
    '''

    betse_cli('sim', betse_sim_config_default.config.filename)
