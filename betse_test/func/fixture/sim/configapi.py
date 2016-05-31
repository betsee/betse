#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes returned by simulation configuration fixtures to other
fixtures and tests requiring a simulation configuration.

These classes provide temporary simulation configurations and metadata on these
configurations for use by public simulation configuration fixtures, which
usually modify the contents of these configurations so as to test specific
feature sets and edge cases.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.wrapper import SimConfigWrapper
from betse.util.path import dirs
from betse.util.type import strs, types
from betse_test.util import requests

# ....................{ CLASSES                            }....................
class SimTestState(object):
    '''
    Simulation configuration context encapsulating simulation configuration,
    state, and metadata.

    Simulation configuration fixtures typically return instances of this class
    as a means of communicating this context to other fixtures and tests.

    Attributes
    ----------
    config : SimConfigWrapper
        Simulation configuration wrapper wrapping the low-level dictionary
        deserialized from the YAML-formatted simulation configuration file with
        path `config_filepath`. Note the contents of this in-memory dictionary
        may be desynchronized from those of this file. For efficiency, callers
        may modify this dictionary to suite test requirements _before_
        reserializing this dictionary back to this file.
    config_filepath : py.path.local
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture as a `py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented `os.path` module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official `py.path` class documentation.
    '''


    def __init__(self, config_filepath: 'py.path.local') -> None:
        '''
        Initialize this simulation configuration context.

        This method copies BETSE's default simulation configuration file,
        complete with all external assets (e.g., geometry masks) referenced and
        required by this file, to the passed path.

        Parameters
        ----------
        config_filepath : py.path.local
            Absolute path to which this method copies BETSE's default simulation
            configuration file. If this file already exists, an exception is
            raised.
        '''
        assert types.is_py_path_local(config_filepath), (
            types.assert_not_py_path_local(config_filepath))

        # Classify the passed parameters. While the "self.config" object
        # classified below provides this filename as a low-level string, this
        # high-level "py.path.local" instance is useful in fixtures and tests.
        self.config_filepath = config_filepath

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = SimConfigWrapper.wrap_new_default(
            filename=str(config_filepath))

        # Unconditionally disable configuration options either requiring
        # interactive input *OR* displaying interactive output for all fixtures
        # and tests.
        self.config.disable_interaction()


    def get_command_context(self) -> 'contextlib.contextmanager':
        '''
        Context manager changing the current working directory (CWD) of the
        current test to the directory containing this configuration file for the
        duration of this context.

        This method is dynamically called by the `run()` method of the
        `CLITester` instance returned by the `betse_cli` fixture for the
        current test.
        '''

        return dirs.current(self.config.dirname)

# ....................{ MAKERS                             }....................
# To force this fixture to return a new object for all parent fixtures, this
# fixture is declared with default scope (i.e., test).

# Ideally, this utility function would be decorated to be a utility fixture --
# permitting parent fixtures calling this function to request this fixture,
# permitting parent fixtures to *NOT* repeatedly request the same fixtures
# required by this function.
#
# Doing so would require this fixture to have function scope. Since parent
# fixtures have session scope, however, py.test prohibits the latter from
# requesting the former. Ergo, this is a utility function instead.
def make(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> 'SimTestState':
    '''
    Create a temporary simulation configuration specific to the parent fixture
    _and_ return a test-specific object encapsulating this configuration.

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
    SimTestState
        Test-specific object encapsulating a temporary simulation configuration
        file specific to the parent fixture.

    See Also
    ----------
    http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
        Official `tmpdir_factory` fixture parameter documentation.
    https://pytest.org/latest/builtin.html#_pytest.python.FixtureRequest
        Official `request` fixture parameter documentation.
    '''

    # Name of the parent fixture.
    parent_fixture_name = requests.get_fixture_name(request)

    # Basename of the temporary directory to be returned for the parent fixture.
    # Equivalently, this is the unique suffix of the name of the parent fixture
    # following the mandatory prefix "betse_sim_config_" in that name. Since the
    # parent fixture does *NOT* physically call this fixture, however, that name
    # cannot be inspected from the call stack. While py.test provides no direct
    # means of obtaining this name, it does provide the "request.fixturenames"
    # list of the names of all fixtures required by the current test, which may
    # then be iteratively searched for the expected name. (See above.)
    sim_config_dir_basename = strs.remove_prefix(
        text=parent_fixture_name,
        prefix='betse_sim_config_',
        exception_message=(
            'Simulation configuration fixture "{}" '
            'not prefixed by "betse_sim_config_".'.format(parent_fixture_name)
        ),
    )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object.
    sim_config_dirpath = tmpdir_factory.mktemp(sim_config_dir_basename)

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Create and return a test-specific object encapsulating this file.
    return SimTestState(config_filepath=sim_config_filepath)
