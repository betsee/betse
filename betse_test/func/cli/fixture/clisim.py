#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixtures.

These fixtures automate testing of the current versions of BETSE's externally
runnable CLI and GUI commands (e.g., `betse`, `betse-qt4`) in subprocesses of
the current test process, regardless of whether these commands have been
editably installed (i.e., as synchronized symlinks rather than desynchronized
copies) into the current Python environment or not.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture
from betse.util.type import strs
from betse.util.type.types import SequenceTypes
from betse_test.func.cli.fixture.cliapi import CLITesterPreArged
from betse_test.func.fixture.sim.configapi import SimTestState
from betse_test.util import requests
from betse_test.util.mark.param import serialize_parametrized_fixture

# ....................{ CONSTANTS                          }....................
_CLI_SIM_SUBCOMMANDS_ARGS = (
    ('seed',),
    ('init',),
    ('sim',),
    ('plot', 'seed'),
    ('plot', 'init'),
    ('plot', 'sim'),
)
'''
List of all argument lists running simulation-specific BETSE CLI subcommands,
used to parametrize the `betse_cli_sim` fixture.

**Order is significant.** These subcommands are run by this fixture in the
listed order. Subcommands requiring the output of prior subcommands as input
must be ordered such that the former follows the latter in this list.
'''


# Dynamically synthesize this list by hyphenating the arguments comprising each
# of the above argument lists.
_CLI_SIM_SUBCOMMANDS_ARGS_IDS = tuple(
    strs.join_on(*args, delimiter='-')
    for args in _CLI_SIM_SUBCOMMANDS_ARGS
)
'''
List of all human-readable unique identifiers to be assigned to each argument
list of the `_CLI_SIM_SUBCOMMANDS_ARGS` global, identifying the parameters
accepted by the `betse_cli_sim` fixture.
'''

# ....................{ FIXTURES                           }....................
# To force these fixtures to return new objects for all parent fixtures and
# tests, these fixtures are declared to have default scope (i.e., test).

@serialize_parametrized_fixture
@fixture(params=_CLI_SIM_SUBCOMMANDS_ARGS, ids=_CLI_SIM_SUBCOMMANDS_ARGS_IDS)
def betse_cli_sim(
    betse_cli: 'CLITester',
    request: '_pytest.python.FixtureRequest',
) -> CLITesterPreArged:
    '''
    Fixture returning an instance of the `CLITestMultiRunner` class, suitable
    for iteratively running _all_ simulation-specific BETSE CLI subcommands
    (e.g., `seed`, `init`, `sim`) with the simulation configuration required by
    the current test or fixture.

    Parameters
    ----------
    betse_cli : CLITester
        Object running a single simulation-specific BETSE CLI subcommand.
    request : _pytest.python.FixtureRequest
        Builtin fixture describing this fixture's parent fixture or test.

    Returns
    ----------
    CLITesterPreArged
        Object running the simulation-specific BETSE CLI subcommand defined by
        the current parametrization with the current simulation configuration.
    '''

    # Name of the simulation configuration fixture required by this test.
    sim_config_fixture_name = requests.get_requested_fixture_name_prefixed_by(
        request=request, fixture_name_prefix='betse_sim_config_')

    # Simulation configuration fixture required by this test.
    sim_state = requests.get_requested_fixture(request, sim_config_fixture_name)
    assert isinstance(sim_state, SimTestState), (
        'Object "{}" not a simulation configuration fixture.'.format(sim_state))

    # Argument list comprising the currently parametrized BETSE CLI subcommand
    # passed the basename of this simulation configuration file, validating that
    # this simulation configuration fixture has changed the current working
    # directory (CWD) to this file's directory.
    subcommand_args = list(requests.get_fixture_param(request))
    subcommand_args.append(sim_state.config.basename)

    # Return a new CLI runner specific to the current test.
    return CLITesterPreArged(
        cli=betse_cli,
        subcommand_args=subcommand_args,
    )
