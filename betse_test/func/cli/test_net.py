#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising all simulation subcommands specific to
biochemical reaction and gene regulatory networks (e.g., `betse sim-brn`, `betse
sim-grn`).
'''

# ....................{ IMPORTS                            }....................
# import pytest
# from betse_test.util.mark.fail import xfail

# ....................{ TESTS                              }....................
def test_cli_brn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default biochemical reaction network (BRN) isolated away
    from all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Test the BRN subcommand and all subcommands required by that subcommand.
    betse_cli_sim.run_subcommands(('seed',), ('sim-brn',),)


def test_cli_grn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default gene regulatory network (GRN) isolated away from
    all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Test the GRN subcommand and all subcommands required by that subcommand.
    betse_cli_sim.run_subcommands(('seed',), ('sim-grn',),)


def test_cli_sim_brn_grn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default biochemical reaction network (BRN) _and_
    gene regulatory network (GRN) alongside all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable these networks.
    betse_cli_sim.sim_state.config.enable_networks()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_default()
