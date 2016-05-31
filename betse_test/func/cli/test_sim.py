#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
# from betse_test.mark.skip import skip

# ....................{ TESTS                              }....................
def test_cli_sim_default(
    betse_cli_sim, betse_sim_config_default) -> None:
    '''
    Test the simulation of the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLITestRunner
        Test-specific object encapsulating the simulation-specific BETSE CLI.
    betse_sim_config : SimTestState
        Test-specific object encapsulating this simulation configuration file.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    betse_cli_sim()
