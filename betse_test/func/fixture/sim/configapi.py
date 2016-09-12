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
from betse.util.type import strs
from betse.util.type.types import type_check, CallableTypes, NoneType
from betse_test.util import requests
from py._path.local import LocalPath

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
    config_filepath : LocalPath
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture as a `py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented `os.path` module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official `py.path` class documentation.
    '''


    @type_check
    def __init__(self, config_filepath: LocalPath) -> None:
        '''
        Initialize this simulation configuration context.

        This method copies BETSE's default simulation configuration file,
        complete with all external assets (e.g., geometry masks) referenced and
        required by this file, to the passed path.

        Parameters
        ----------
        config_filepath : LocalPath
            Absolute path to which this method copies BETSE's default simulation
            configuration file as a `py.path.local` instance. If this file
            already exists, an exception is raised.
        '''

        # Classify the passed parameters. While the "self.config" object
        # classified below provides this filename as a low-level string, this
        # high-level "py.path.local" instance is useful in fixtures and tests.
        self.config_filepath = config_filepath

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = SimConfigWrapper.wrap_new_default(
            filename=str(config_filepath))

        # For all child fixtures and tests, unconditionally:
        #
        # * Disable configuration options either requiring interactive input
        #   *OR* displaying interactive output.
        # * Minimize the space and time costs associated with running the
        #   simulation configured by this configuration while preserving all
        #   fundamental configuration features.
        self.config.disable_interaction()
        self.config.minify()


    def get_command_context(self) -> 'contextlib.contextmanager':
        '''
        Context manager changing the current working directory (CWD) of the
        current test to the directory containing this configuration file for the
        duration of this context.

        This method is dynamically called by the `run()` method of the
        `CLITester` instance provided by the `betse_cli` fixture for the current
        fixture or test.
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
#
# Test suite insanity prevails.
@type_check
def make(
    request, tmpdir_factory,
    sim_config_dir_basename: (str, NoneType,) = None,
    config_modifier: CallableTypes + (NoneType,) = None,
) -> SimTestState:
    '''
    Create a temporary simulation configuration specific to the parent fixture,
    call the passed callable accepting (and presumably modifying) this
    configuration if non-`None`, and return a test object encapsulating this
    configuration.

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
    sim_config_dir_basename : optional[str]
        Basename of the temporary directory containing this simulation
        configuration to be created. Defaults to `None`, in which case this
        basename defaults to the name of the current test excluding the prefix
        `test_` (e.g., `cli_sim_default` for the `test_cli_sim_default` test).
    config_modifier : optional[CallableType]
        Callable (e.g., function, lambda, method) accepting the newly created
        simulation configuration dictionary as a mandatory parameter. For
        convenience, this configuration is implicitly written back to disk
        _after_ this callable is called. This callable is merely a shorthand
        convenience for fixtures modifying this configuration in a trivial
        manner expressible as a `lambda` statement. Fixtures requiring less
        trivial modifications may do so by directly modifying and rewriting the
        `config` attribute of the returned object instead. Defaults to `None`,
        in which case this parameter is ignored.

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

    # If no basename was passed...
    if sim_config_dir_basename is None:
        # Name of the current test.
        test_name = requests.get_tested_name(request)
        # print('    request.node: {}'.format(request.node))
        # print('    test_name: {}'.format(test_name))

        # Default this basename to the name of the current test excluding the
        # prefix "test_".
        sim_config_dir_basename = strs.remove_prefix(
            text=test_name,
            prefix='test_',
            exception_message='Test name "{}" not prefixed by "test_".'.format(
                test_name),
        )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object.
    sim_config_dirpath = tmpdir_factory.mktemp(sim_config_dir_basename)

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Test-specific object encapsulating this simulation configuration file.
    sim_state = SimTestState(config_filepath=sim_config_filepath)

    # If the calling fixture requested a configuration modification, do so.
    if config_modifier is not None:
        config_modifier(sim_state.config)

    # Unconditionally write this configuration to disk *AFTER* possibly
    # modifying this configuration.
    sim_state.config.overwrite()

    # Return this encapsulation.
    return sim_state
