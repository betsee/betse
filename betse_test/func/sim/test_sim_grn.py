#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising all simulation subcommands pertaining
to gene regulatory networks (e.g., `betse sim-grn`, `betse plot sim-grn`).
'''

# ....................{ IMPORTS                            }....................
import pytest
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
@pytest.mark.parametrize(
    ('unpickle_phase_name', 'init_subcommands'), (
        ('seed', (('seed',),),),
        ('init', (('seed',), ('init',),),),
        ('sim',  (('seed',), ('init',), ('sim',),),),
    ))
def test_cli_grn_isolated(
    betse_cli_sim: 'CLISimTester',
    unpickle_phase_name: str,
    init_subcommands: tuple,
) -> None:
    '''
    Test simulating the default gene regulatory network (GRN) isolated away from
    all bioelectrical phenomena for the passed unpickle simulation phase and
    simulation subcommands.

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    unpickle_phase_name : str
        Lowercase name of the type of **unpickle simulation phase** (i.e.,
        previously pickled simulation phase to unpickle as the computational
        basis for the current network to be run by the ``betse sim-grn``
        subcommand) to be exercised by this test.
    init_subcommands : tuple
        Tuple of all simulation subcommands to be exercised by this test
        *before* exercising the ``betse sim-grn`` subcommand.
    '''

    # Defer heavyweight imports required for logging. For readability in the
    # unfortunate event of exceptions raised by subsequent imports, defer all
    # remaining imports.
    from betse.util.io.log import logs

    # Log this networking attempt with suitable aesthetics.
    logs.log_banner(
        title='sim-grn ({})'.format(unpickle_phase_name), padding='~')

    # Defer all remaining heavyweight imports.
    from betse.science.config.confenum import GrnUnpicklePhaseType
    from betse.util.path import pathnames
    from betse.util.type import enums

    # Simulation configuration specific to this test.
    p = betse_cli_sim.sim_state.p

    # Enable the saving of visuals, preventing the "plot sim-grn" subcommand
    # tested below from silently reducing to a noop.
    betse_cli_sim.sim_state.config.enable_visuals_save()

    #FIXME: Replace this overkill method call with direct usage of the local
    #"p" parameter; then remove this method entirely from the codebase.
    # Enable networking.
    betse_cli_sim.sim_state.config.enable_networks()

    # Enable this type of unpickle simulation phase by mapping from the passed
    # string to the corresponding enumeration member. Ideally, this member
    # rather than this string would simply be passed to this test;
    # unfortunately, doing so would require unsafely importing from the main
    # codebase at module scope, encouraging non-human-readable exceptions.
    p.grn_unpickle_phase_type = enums.get_member_from_name_uppercased(
        enum_type=GrnUnpicklePhaseType, enum_member_name=unpickle_phase_name)

    # Test all simulation subcommands required by this GRN-specific subcommand
    # with this configuration for this phase type.
    betse_cli_sim.run_subcommands(*init_subcommands)

    # Test this GRN-specific subcommand from scratch.
    betse_cli_sim.run_subcommands(('sim-grn',),)

    # Log this re-networking attempt with suitable aesthetics.
    logs.log_banner(
        title='sim-grn (re-{})'.format(unpickle_phase_name), padding='~')

    # Prepare to rerun the "sim-grn" subcommand from the prior run pickled by
    # the prior subcommand. To do so safely (in order):
    #
    # * Reconstruct the relative filename of the prior GRN run.
    # * Ensure that the next GRN run is pickled to another file.
    p.grn_unpickle_filename_relative = pathnames.join(
        p.grn_pickle_dirname_relative, p.grn_pickle_basename)
    p.grn_pickle_basename = 'new_' + p.grn_pickle_basename

    # Redefine all absolute pathnames depending upon these relative pathnames.
    p.reload_paths()

    # Test both this GRN-specific subcommand from the prior such run and, for
    # completeness, exporting the results of doing so.
    betse_cli_sim.run_subcommands(('sim-grn',), ('plot', 'sim-grn',),)


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
