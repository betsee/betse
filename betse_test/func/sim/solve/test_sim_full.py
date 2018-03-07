#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising the full simulation solver.
'''

# ....................{ TESTS                              }....................
def test_cli_sim_full_noecm(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Functional test exporting all available exports (e.g., CSVs, plots,
    animations) with all simulation features required by these exports,
    including the full solver but excluding extracellular spaces.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable all exports and features required by these exports, excluding ECM.
    betse_cli_sim.sim_state.config.enable_solver_full_exports_noecm()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_try()


def test_cli_sim_full_ecm(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Functional test exporting all available exports (e.g., CSVs, plots,
    animations) with all simulation features required by these exports,
    including both the full solver and extracellular spaces.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable all exports and features required by these exports, including ECM.
    betse_cli_sim.sim_state.config.enable_solver_full_exports_ecm()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_try()


def test_cli_sim_full_vg_ions(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Functional test simulating all voltage-gated ion channels (e.g., sodium,
    potassium) *and* simulation features required by these channels, including
    the full solver.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable all voltage-gated ion channels and features required by these
    # channels.
    betse_cli_sim.sim_state.config.enable_solver_full_vg_ions()

    # Test all simulation-specific subcommands *EXCLUDING* plotting subcommands
    # (which other tests already exercise) with this configuration.
    betse_cli_sim.run_subcommands_sim()
