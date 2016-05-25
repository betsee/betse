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

#FIXME: We've added a new "--dist-dir" option to our "freeze_*" family of
#setuptools subcommands. Now, let's add support for fixtures running these
#subcommands. Is there any means of programmatically running setuptools
#subcommands in the current Python process? We suspect not, sadly.

#FIXME: Add seemless support for exercising all CLI-specific functional tests
#leveraging the "betse_cli" fixture against a frozen rather than unfrozen
#version of BETSE's CLI. Do *NOT* parametrize this fixture to unconditionally
#run each such test against both frozen and unfrozen versions of BETSE's CLI, as
#doing so would substantially increase space and time complexity with little to
#no tangible gain. Instead:
#
#* Define a new "--is-frozen" command-line option specific to our setuptools
#  "test" subcommand, defaulting to disabled.
#* When this option is *NOT* passed:
#  * One and only one PyInstaller test should be run -- say, test_pyi_try().
#    This test is intended only to provide a coarse-grained sanity check of
#    whether or not:
#    * BETSE is actually freezable under PyInstaller.
#    * The resulting executable is successfully simulatable with the default
#      simulation configuration via "betse try".
#  * Consequently, this test should (in order):
#    1. Freeze BETSE in one-dir mode (which requires no additional compression
#       and decompression steps and hence is slightly more time efficient for
#       our simple purposes) into a temporary subdirectory.
#    2. Run "betse try" against this frozen executable.
#  * All other CLI tests should be run as is against the unfrozen version of
#    BETSE importable by the active Python interpreter.
#* When this option is passed:
#  * The aforementioned test_pyi_try() test should *NOT* be run, as doing so
#    would be wholly redundant.
#  * BETSE should be frozen in one-dir mode, as discussed above.
#  * All other CLI tests should be run against this frozen executable rather
#    than the unfrozen, importable version of BETSE.
#
#Hence, the "--is-frozen" CLI option serves as a high-level switch fundamentally
#changing testing behaviour. The intention, of course, is that this option would
#only be passed immediately before producing an official new frozen version of
#BETSE. This provides a sanity check on frozen executable behaviour otherwise
#difficult (if not infeasible) to do manually.
#FIXME: This is great. We'd only like to make one improvement to the
#aforementioned PyInstaller testing support: marks. Use them. Explicit is
#typically better than implicit. In particular:
#
#* Define a new "betse_test.mark.able" submodule.
#* Define a new "@freezable" mark decorator in this submodule.
#* Decorate all tests suitable for running under a frozen copy of BETSE with
#  this decorator. This should probably *ONLY* include:
#  * All functional CLI tests, excluding "test_pyi_try". Everything else (e.g.,
#    unit tests) are *NOT* suitable for running frozen.
#FIXME: Urgh. That's great, but it's still not quite enough. Why? Simulation
#configurations. We currently create them by importing BETSE packages, which
#clearly won't work for frozen commands. Instead, if running frozen, we'll need
#to just call "betse config". Simple, if tedious.

# ....................{ IMPORTS                            }....................
from betse.util.type import objects
from betse_test.func.cli.fixture.commandapi import CLITestRunner
from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
# To force these fixtures to return new objects for all parent fixtures and
# tests, these fixtures is declared to have default scope (i.e., test).

#FIXME: Implement me. See commentary preceding test_cli_sim_default().

@fixture
def betse_cli_sim(request: '_pytest.python.FixtureRequest') -> CLITestRunner:
    '''
    Fixture returning an instance of the `CLITestRunner` class, suitable for
    iteratively running _all_ simulation-specific BETSE CLI subcommands (e.g.,
    `seed`, `init`) under the simulation configuration required by the current
    fixture or test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).

    Returns
    ----------
    CLITestRunner
        Object running _all_ simulation-specific BETSE CLI subcommands.
    '''

    raise ValueError('Not implemented yet.')


@fixture
def betse_cli(request: '_pytest.python.FixtureRequest') -> CLITestRunner:
    '''
    Fixture returning an instance of the `CLITestRunner` class, suitable for
    running the BETSE CLI command required by the current fixture or test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).

    Returns
    ----------
    CLITestRunner
        Object running the BETSE CLI command.
    '''

    # Names of all BETSE-specific fixtures required by the current test,
    # excluding the current fixture.
    betse_fixture_names = requests.get_fixture_names_prefixed_by(
        request=request, fixture_name_prefix='betse_')

    # List of the command-specific context managers provided by these fixtures,
    # corresponding to the context managers returned by the optional
    # get_command_context() methods defined by the instances of these fixtures.
    betse_fixture_command_contexts = []

    # For the name of each such fixture...
    for betse_fixture_name in betse_fixture_names:
        # Fixture object returned by this fixture.
        betse_fixture = requests.get_fixture(
            request=request, fixture_name=betse_fixture_name)

        # get_command_context() method defined by this fixture object if any.
        betse_fixture_get_command_context = objects.get_method_or_none(
            obj=betse_fixture, method_name='get_command_context')

        # If this fixture object defines this method, append the context manager
        # returned by call this method to this list.
        if betse_fixture_get_command_context is not None:
            betse_fixture_command_contexts.append(
                betse_fixture_get_command_context())

    # Return a new CLI runner specific to the current test.
    return CLITestRunner(contexts=betse_fixture_command_contexts)
