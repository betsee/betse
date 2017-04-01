#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising all simulation subcommands specific to
biochemical reaction and gene regulatory networks (e.g., `betse sim-brn`, `betse
sim-grn`).
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.skip import skip_unless_lib_runtime_optional

# ....................{ GLOBALS                            }....................
# Decorator skipping all tests running network plotting subcommands (e.g., "plot
# sim-brn", "plot sim-grn") requiring these optional runtime dependencies.
skip_unless_networkable = skip_unless_lib_runtime_optional('networkx', 'pydot')

# ....................{ TESTS                              }....................
#FIXME: Enable plots and animations for these tests.
@skip_unless_networkable
def test_cli_brn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default biochemical reaction network (BRN) isolated away
    from all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable the saving of visuals, preventing the "plot sim-brn" subcommand
    # tested below from silently reducing to a noop.
    betse_cli_sim.sim_state.config.enable_visuals_save()

    # Test all BRN-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands(
        ('seed',), ('sim-brn',), ('plot', 'sim-brn',),)


@skip_unless_networkable
def test_cli_grn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default gene regulatory network (GRN) isolated away from
    all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable the saving of visuals, preventing the "plot sim-grn" subcommand
    # tested below from silently reducing to a noop.
    betse_cli_sim.sim_state.config.enable_visuals_save()

    # Test all GRN-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands(
        ('seed',), ('sim-grn',), ('plot', 'sim-grn',),)


@skip_unless_networkable
def test_cli_sim_brn_grn(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default biochemical reaction network (BRN) *and*
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
