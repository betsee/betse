#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
from betse_test.mark.skip import skip

# ....................{ TESTS                              }....................
#FIXME: While useful, this should be split apart into its constituent steps,
#most of which should probably become fixtures for efficient reuse. Wait, no.
#Ideally, we should construct a linear chain of tests, each exercising a
#simulation phase (e.g., "seed", "init", "sim") depending on the success of the
#prior test in this chain.
#FIXME: Right. No. This is clearly a job for parametrization. Ideally, we only
#want to declare a single test for each simulation configuration which accepts
#the names of all CLI subcommands to be passed the configuration file. Perhaps
#we want to replace the "betse_cli" fixture requested below with a new
#"betse_cli_sim" fixture producing a class inheriting the class produced by the
#former and adding support for parametrization? Consider it.
#
#Indeed, that's the way. However, we don't want tests themselves to specify
#parameters. Instead, we want the new "betse_cli_sim" fixture to unconditionally
#parametrize *ITSELF* by the in-order list of all CLI subcommands to be run:
#e.g.,
#
#    (
#        ('seed'), ('init'), ('sim'),
#        ('plot', 'seed'), ('plot', 'init'), ('plot', 'sim'),
#    )
#
#Why? A couple obvious reasons:
#
#1. We don't want to duplicate that list for every test. By centralizing that
#   list into a single fixture, changes are dramatically simplified.
#2. Tests may *THINK* they know which CLI subcommands they need run, but...
#   really, they don't. Each new configuration really requires everything to be
#   revalidated from scratch. This approach forces that.
#3. By parametrizing the "betse_cli_sim" fixture, all tests requesting that
#   fixture will themselves receive the same parameters -- which means that
#   users can interactively test exact CLI subcommands for exact configurations
#   by simply using "./test -k". Great!
#
#In short, this is the way forward.

def test_cli_sim_default(
    betse_cli, betse_sim_config_default) -> None:
    '''
    Test the simulation of the default simulation configuration.

    Parameters
    ----------
    betse_cli : CLITestRunner
        Test-specific object encapsulating the BETSE CLI.
    betse_sim_config : SimTestConfig
        Test-specific object encapsulating this simulation configuration file.
    '''

    # Pass the basename rather than absolute path of this configuration file,
    # testing that the "betse_sim_config_default" fixture has set the current
    # working directory (CWD) to the directory containing this file.
    betse_cli('sim', betse_sim_config_default.config.basename)
