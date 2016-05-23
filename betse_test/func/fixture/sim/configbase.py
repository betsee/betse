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
from betse_test.util import requests
from betse_test.util.exceptions import BetseTestFixtureException
from pytest import fixture

# ....................{ FIXTURES                           }....................
#FIXME: Conditionally delete all temporary directories if and only if all
#tests requiring those directories succeeded. Since we have idea how to safely
#implement this, at the moment, temporary directories are *NOT* currently
#automatically deleted. When the are, add the following documentation to the
#fixture docstring below:
#
#    "This temporary directory will be recursively deleted on completion of this
#     test session _without_ requiring manual intervention (e.g., via
#     finalizers) from any other fixtures or tests."
#
#The simplest means of implementing this would probably be to redeclare the
#"_betse_with_sim_config" fixture using the @pytest.yield_fixture decorator,
#which is bar-none the simplest means of implementing this sort of "teardown"
#logic safely and efficiently. See:
#
#    http://pytest.org/latest/yieldfixture.html

# To force this fixture to return a new object for all parent fixtures, this
# fixture is declared with default scope (i.e., test).
@fixture
def _betse_with_sim_config(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> 'SimTestConfig':
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

    This fixture temporarily changes the current directory of the active Python
    process to this temporary directory. Default simulation configuration paths
    are relative to the directory containing the simulation configuration file:
    namely, this temporary directory. Changing directories resolves these paths
    to this directory. (Failing to do so would incorrectly resolve these paths
    to the current directory, with predictably disastrous outcomes.) While this
    fixture could instead globally search-and-replace all relative simulation
    configuration paths with absolute paths, doing so would be considerably more
    complex, fragile, and error-prone than simply changing directories.

    This temporary directory and hence simulation configuration will be
    accessible for the duration of this test session. Tests and fixtures may
    reuse this configuration, but should do so _only_ in an orderly and
    coordinated manner preserving deterministic testing.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin session-scoped fixture whose `mktemp()` method returns a
        `py.path.local` instance encapsulating a new temporary directory.

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
    # prefixed by "betse_sim_config_". Since this is a session-scope fixture
    # applicable to multiple tests, inspecting the "requests" object specific to
    # a single test is usually unsafe. In this case, however, the following
    # logic only finds the parent fixture commonly required by all tests
    # requiring that fixture including the first such test, and hence is safe.
    #
    # See below for additional discussion.
    sim_config_fixture_names = [
        fixture_name
        for fixture_name in requests.get_test_fixture_names(request)
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
