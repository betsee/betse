#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures and fixture classes creating temporary simulation configurations
isolated to specific tests, which typically modify the contents of these
configurations so as to exercise specific feature sets and edge cases.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.wrapper import SimConfigWrapper
from betse.util.path import dirs
from betse.util.type import strs
from betse.util.type.types import type_check
from betse_test.util import requests
from py._path.local import LocalPath
from pytest import fixture

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


    def context(self) -> 'contextlib.contextmanager':
        '''
        Context manager changing the current working directory (CWD) of the
        current test to the directory containing this configuration file for the
        duration of this context.

        Default simulation configuration paths are relative to the directory
        containing the simulation configuration file: namely, this temporary
        directory. Changing directories resolves these paths to this directory.
        (Failing to do so would incorrectly resolve these paths to the current
        directory, with predictably disastrous outcomes.) While this class
        could instead globally search-and-replace all relative simulation
        configuration paths with absolute paths, doing so would be considerably
        more complex, fragile, and error-prone than simply changing directories.
        '''

        # Defer to the generator returned by the following utility function.
        return dirs.current(self.config.dirname)

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_sim_config(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> SimTestState:
    '''
    Per-test fixture creating a temporary default simulation configuration file
    and returning an object encapsulating the contents of this file.

    Configuration Modifications (On-disk)
    ----------
    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory whose basename is the name
    of the test requesting this fixture excluding the prefixing substring
    `test_`. As example, when requested by the `test_cli_sim_default` test,
    this fixture creates a temporary simulation configuration file
    `{tmpdir}/cli_sim_default/sim_config.yaml` (e.g.,
    `/tmp/pytest-0/cli_sim_default/sim_config.yaml`), where `{tmpdir}` is the
    absolute path of this test session's root temporary directory (e.g.,
    `/tmp/pytest-0/`).

    This temporary directory and hence simulation configuration will be
    accessible _only_ for the duration of the current test. Subsequently run
    tests and fixtures may _not_ safely reuse this configuration.

    Configuration Modifications (In-memory)
    ----------
    This fixture also transforms the in-memory instance of the
    :class:`betse.science.config.wrapper.SimConfigWrapper` class encapsulating
    this configuration as follows:

    * All configuration options either requiring interactive input _or_
      displaying interactive output are disabled (e.g., plots, animations).
    * The space and time costs associated with simulating this configuration
      are safely minimized in a manner preserving all features.

    Since this fixture does _not_ write these changes back to this file, the
    parent fixture or test is expected to do so manually (e.g., by calling the
    :meth:`betse.science.config.wrapper.SimConfigWrapper.overwrite` method on
    the `config` attribute of the object returned by this fixture).

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin session-scoped fixture whose `mktemp()` method returns a
        :class:`py.path.local` instance encapsulating a new temporary directory.

    Returns
    ----------
    SimTestState
        Test-specific object encapsulating a temporary simulation configuration
        file specific to the current test, including such metadata as:
        * The absolute path of this configuration's on-disk YAML file.
        * This configuration's in-memory dictionary deserialized from this file.
    '''

    # Name of the current test.
    test_name = requests.get_tested_name(request)
    # print('    request.node: {}'.format(request.node))
    # print('    test_name: {}'.format(test_name))

    # Basename of the temporary directory containing this configuration file,
    # set to the name of the current test excluding the prefixing "test_".
    sim_config_dir_basename = strs.remove_prefix(
        text=test_name,
        prefix='test_',
        exception_message=(
            'Test name "{}" not prefixed by "test_".'.format(test_name)),
    )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object. See also:
    #     http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
    sim_config_dirpath = tmpdir_factory.mktemp(sim_config_dir_basename)

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Test-specific object encapsulating this simulation configuration file.
    sim_state = SimTestState(config_filepath=sim_config_filepath)

    # Return this object.
    return sim_state
