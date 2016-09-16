#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_cli_sim_default(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running multiple BETSE CLI simulation subcommands.
    '''

    # Persist all in-memory configuration changes back to disk.
    betse_cli_sim.sim_state.config.overwrite()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_default()


def test_cli_sim_visuals(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating all exported visuals (e.g., plots, animations) _and_
    simulation features required by these visuals.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running multiple BETSE CLI simulation subcommands.
    '''

    # Enable all exported visuals and features required by these visuals.
    betse_cli_sim.sim_state.config.enable_visuals()

    # Persist all in-memory configuration changes back to disk.
    betse_cli_sim.sim_state.config.overwrite()
    # print('\n!!!!!after solving: {}'.format(betse_cli_sim.sim_state.config._config['results options']['after solving']))

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_default()


#FIXME: Convert into a functional test of the above format.
# @fixture(scope='session')
# def betse_sim_config_video(
#     request: '_pytest.python.FixtureRequest',
#     tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
# ) -> SimTestState:
#     '''
#     Context manager-driven fixture creating a temporary simulation configuration
#     enabling encoding of at least one animation and all features required by
#     this animation as compressed video of the passed filetype with the preferred
#     matplotlib animation writer of the passed name _and_ returning a
#     test-specific object encapsulating this configuration.
#
#     See Also
#     ----------
#     :func:`betse_sim_config_anims`
#         Further details on fixture parameters and return value.
#     '''
#
#     return configapi.make(
#         request, tmpdir_factory,
#         config_modifier=lambda config: config.enable_video())
