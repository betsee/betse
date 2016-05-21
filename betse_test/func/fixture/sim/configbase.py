#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Private simulation configuration fixtures.

These fixtures provide temporary simulation configurations and metadata on these
configurations for use by public simulation configuration fixtures, which
usually modify the contents of these configurations so as to test specific
feature sets and edge cases.

These fixtures are _not_ intended to be directly required by tests.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture
from betse_test.exceptions import BetseTestFixtureException

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def _betse_sim_config(tmpdir_factory, request) -> 'SimTestConfig':
    '''
    Context manager-driven fixture creating a temporary simulation configuration
    specific to the parent fixture _and_ returning a test-specific object
    encapsulating this configuration.

    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory whose basename is the
    name of the parent fixture excluding the prefixng substring
    `betse_sim_config_`. If the name of the parent fixture is _not_ prefixed by
    `betse_sim_config_` or this fixture is _not_ called by another fixture, an
    exception is raised. For example, when required by the public parent fixture
    `betse_sim_config_default`, this fixture creates a temporary simulation
    configuration file `{tmpdir}/default/sim_config.yaml` (e.g.,
    `/tmp/pytest-0/default/sim_config.yaml`), where `{tmpdir}` is the absolute
    path of this session's root temporary directory (e.g., `/tmp/pytest-0/`).

    This temporary directory and hence simulation configuration will be
    accessible for the duration of this test session. Tests and fixtures may
    reuse this configuration, but should do so _only_ in an orderly and
    coordinated manner preserving deterministic testing.

    This temporary directory will be recursively deleted on completion of this
    test session _without_ requiring manual intervention (e.g., via finalizers)
    from any other fixtures or tests.

    Parameters
    ----------
    tmpdir_factory : ???
        Builtin session-scoped fixture whose `mktemp()` method returns a
        `py.path.local` instance encapsulating a new temporary directory.
    request : FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).

    Returns
    ----------
    SimTestConfig
        Test-specific object encapsulating a temporary simulation configuration
        file specific to the parent fixture.

    See Also
    ----------
    http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
        Official `tmpdir_factory` fixture parameter documentation.
    https://pytest.org/latest/builtin.html#_pytest.python.FixtureRequest
        Official `request` fixture parameter documentation.
    '''

    # Defer imports importing heavyweight dependencies to their point of use.
    from betse_test.func.fixture.sim.configapi import SimTestConfig

    # Names of all fixtures required by the current test whose names are
    # prefixed by "betse_sim_config_". (See below for comments.)
    sim_config_fixture_names = [
        fixture_name
        for fixture_name in request.fixturenames
         if fixture_name.startswith('betse_sim_config_')
    ]

    # Number of such fixtures.
    sim_config_fixture_count = len(sim_config_fixture_names)

    # If this fixture is required by either more or less than exactly one such
    # fixture, raise an exception.
    if sim_config_fixture_count != 1:
        # Exception message to be raised.
        exception_message = None
        if sim_config_fixture_count == 0:
            exception_message = (
                '"_betse_sim_config" fixture not required by '
                '"betse_sim_config_"-prefixed fixture.'
            )
        else:
            exception_message = (
                '"_betse_sim_config" fixture required by multiple '
                '"betse_sim_config_"-prefixed fixtures {}.'.format(
                    sim_config_fixture_names)
            )

        # Raise this exception with this message.
        raise BetseTestFixtureException(exception_message)

    # Basename of the temporary directory to be returned for the parent fixture.
    # Equivalently, this is the unique suffix of the name of the parent fixture
    # following the mandatory prefix "betse_sim_config" in that name. Since the
    # parent fixture does *NOT* physically call this fixture, however, that name
    # cannot be inspected from the call stack. While py.test provides no direct
    # means of obtaining this name, it does provide the "request.fixturenames"
    # list of the names of all fixtures required by the current test, which may
    # then be iteratively searched for the expected name. (See above.)
    sim_config_dir_basename = sim_config_fixture_names[0]

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object.
    sim_config_dirpath = tmpdir_factory.mktemp(sim_config_dir_basename)

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Create and return a test-specific object encapsulating this file.
    return SimTestConfig(config_filepath=sim_config_filepath)
