#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising the equivalent circuit-based simulation
solver.
'''

# ....................{ IMPORTS                            }....................
# from betse_test.util.mark.fail import xfail

# ....................{ TESTS                              }....................
def test_cli_sim_circuit(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Functional test exporting all available exports (e.g., CSVs, plots,
    animations) with all simulation features required by these exports,
    including the equivalent circuit-based solver but excluding extracellular
    spaces (which are only partially supported by this solver).

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable all exports and features required by these exports.
    betse_cli_sim.sim_state.config.enable_solver_circuit_exports()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_try()
