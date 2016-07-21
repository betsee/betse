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
    betse_cli_sim: 'CLITester',
    betse_sim_config_default: 'SimTestState',
) -> None:
    '''
    Test simulating with the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLITester
        Object encapsulating CLI-driven simulation testing.
    betse_sim_config_default : SimTestState
        Object encapsulating the default simulation configuration.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    betse_cli_sim()


def test_cli_sim_anims(
    betse_cli_sim: 'CLITester',
    betse_sim_config_anims: 'SimTestState',
) -> None:
    '''
    Test simulating all animations and features required by these animations.

    Parameters
    ----------
    betse_cli_sim : CLITester
        Object encapsulating CLI-driven simulation testing.
    betse_sim_config_anims : SimTestState
        Object encapsulating the simulation configuration enabling animations.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    betse_cli_sim()
