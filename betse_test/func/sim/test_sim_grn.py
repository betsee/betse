#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising all simulation subcommands pertaining
to gene regulatory networks (e.g., `betse sim-grn`, `betse plot sim-grn`).
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.skip import skip_unless_lib_runtime_optional

# ....................{ DECORATORS                         }....................
skip_unless_networkable = skip_unless_lib_runtime_optional('networkx', 'pydot')
'''
Decorator skipping the decorated test if either of these optional runtime
dependencies are unavailable, both of which are required by network plotting
subcommands (e.g., ``plot sim-grn``).
'''

# ....................{ TESTS                              }....................
@skip_unless_networkable
def test_cli_grn_isolated(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default gene regulatory network (GRN) isolated away from
    all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable these networks.
    betse_cli_sim.sim_state.config.enable_networks()

    # Enable the saving of visuals, preventing the "plot sim-grn" subcommand
    # tested below from silently reducing to a noop.
    betse_cli_sim.sim_state.config.enable_visuals_save()

    # Test all GRN-specific subcommands required to simulate from scratch with
    # this configuration.
    betse_cli_sim.run_subcommands(('seed',), ('sim-grn',),)

    # Test all GRN-specific subcommands required to simulate by reusing the GRN
    # simulated by previous subcommands with this configuration.

    #FIXME: Uncomment this line and comment the subsequent line. Doing so will
    #require setting the "betse_cli_sim.sim_state.p.grn_loadfrom" variable to
    #the concatenation of "betse_cli_sim.sim_state.p.grn_savedir" +
    #"betse_cli_sim.sim_state.p.grn_savefile". Naturally, doing that will in
    #turn necessitate refactoring these variables into YAML data descriptors.

    # betse_cli_sim.run_subcommands(('sim-grn',), ('plot', 'sim-grn',),)
    betse_cli_sim.run_subcommands(('plot', 'sim-grn',),)


@skip_unless_networkable
def test_cli_sim_grn_integrated(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Test simulating the default gene regulatory network (GRN) integrated
    together with all bioelectrical phenomena.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Enable these networks.
    betse_cli_sim.sim_state.config.enable_networks()

    # Test all default simulation-specific subcommands with this configuration.
    betse_cli_sim.run_subcommands_try()
