#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures and fixture classes creating temporary simulation configurations
isolated to specific tests, which typically modify the contents of these
configurations so as to exercise specific feature sets and edge cases.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check
from pytest import fixture
from py._path.local import LocalPath

# ....................{ CLASSES                            }....................
#FIXME: Most use of the increasingly obsolete "SimTestState.config" wrapper
#attribute (both here and everywhere else) should be replaced by use of the new
#"SimTestState.p" property, which increasingly provides all test functionality.

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
        path :attr:`conf_filepath`. Note the contents of this dictionary may
        be desynchronized from those of this file. For efficiency, callers may
        modify this dictionary to suite test requirements *before* reserializing
        this dictionary back to this file.

    Attributes (Path)
    ----------
    conf_dirname : str
        Absolute path of the parent directory containing a temporary simulation
        configuration file specific to the parent fixture, equivalent to
        ``pathnames.get_dirname(str(self.conf_filepath))``.
    conf_filename : str
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture, equivalent to ``str(self.conf_filepath)``.
    conf_filepath : LocalPath
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture as a :class:`py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented :mod:`os.path`
        module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official :mod:`py.path` submodule documentation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, conf_filepath: LocalPath) -> None:
        '''
        Initialize this simulation configuration context.

        This method (in order):

        #. Copies BETSE's default simulation configuration file, complete with
           all external assets (e.g., geometry masks) referenced and required by
           this file, to the passed path.
        #. Sanitizes the copied simulation configuration file for all child
           fixtures and tests by unconditionally:
           * Disabling configuration options either:
             * Requiring interactive input.
             * Displaying interactive output.
           * Minimizing the space and time costs associated with running
             simulations configured by this configuration while preserving all
             fundamental configuration features.

        Caveats
        ----------
        For efficiency, note that the configuration changes applied by this
        method (listed above) reside *only* in-memory; they have yet to be
        written back to disk. Callers are required to do so manually (e.g., the
        ``is_overwriting_config`` parameter passed to the
        :meth:`CLISimTester.run_subcommands` method).

        Parameters
        ----------
        conf_filepath : LocalPath
            Absolute path to which this method will copy the default simulation
            configuration file as a :class:`py.path.local` instance. If this
            file already exists, an exception is raised.
        '''

        # Defer heavyweight imports. This subclass inherits a class defined by
        # the main codebase and is hence *NOT* safely importable above.
        from betse.util.path import pathnames
        from betse_test.fixture.simconfig.simconfwrapper import (
            SimConfigTestWrapper)

        # Classify the passed parameters. While the "self.config" object
        # classified below provides this filename as a low-level string, this
        # high-level "py.path.local" instance is useful in fixtures and tests.
        self.conf_filepath = conf_filepath
        self.conf_filename = str(conf_filepath)
        self.conf_dirname = pathnames.get_dirname(self.conf_filename)

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = SimConfigTestWrapper.wrap_new_default(
            filename=self.conf_filename)

        # Sanitize this configuration for all child fixtures and tests.
        self.config.disable_interaction()
        self.config.minify()

    # ..................{ CONTEXTS                           }..................
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

        # Defer heavyweight imports.
        from betse.util.os.shell import shelldir

        # Defer to the generator returned by the following utility function.
        return shelldir.setting_cwd(self.config.dirname)

    # ..................{ PROPERTIES                         }..................
    # For safety, these properties lack setters and hence are read-only.

    @property
    def p(self) -> 'Parameters':
        '''
        High-level simulation configuration encapsulated by this test wrapper.
        '''

        return self.config.p

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_sim_config(betse_temp_dir: LocalPath) -> SimTestState:
    '''
    Per-test fixture creating a temporary default simulation configuration file
    and returning an object encapsulating the contents of this file.

    Configuration Modifications (On-disk)
    ----------
    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory whose basename is the name
    of the test requesting this fixture excluding the prefixing substring
    ``test_``. When requested by the ``test_cli_sim_default`` test, for example,
    this fixture creates a temporary simulation configuration file
    ``{tmpdir}/cli_sim_default/sim_config.yaml`` for the absolute path
    ``{tmpdir}`` of this test session's root temporary directory (e.g.,
    ``/tmp/pytest-0/cli_sim_default/sim_config.yaml``).

    This directory and hence simulation configuration is safely accessible
    *only* for the duration of the current test. Subsequently run tests and
    fixtures *cannot* safely reuse this configuration.

    Configuration Modifications (In-memory)
    ----------
    This fixture also transforms the in-memory instance of the
    :class:`betse.science.parameters.Parameters` class encapsulating this
    configuration as follows:

    * All configuration options either requiring interactive input *or*
      displaying interactive output are disabled (e.g., plots, animations).
    * The space and time costs associated with simulating this configuration
      are safely minimized in a manner preserving all features.

    Since this fixture does _not_ write these changes back to this file, the
    parent fixture or test is expected to do so manually (e.g., by calling the
    :meth:`SimTestState.config.overwrite` method on the object returned by this
    fixture).

    Parameters
    ----------
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current test.

    Returns
    ----------
    SimTestState
        Test-specific object encapsulating a temporary simulation configuration
        file specific to the current test, including such metadata as:
        * The absolute path of this configuration's on-disk YAML file.
        * This configuration's in-memory dictionary deserialized from this file.
    '''

    # Absolute path of this configuration file in this temporary directory.
    sim_conf_filepath = betse_temp_dir.join('sim_config.yaml')

    # Test-specific object encapsulating this simulation configuration file.
    sim_state = SimTestState(conf_filepath=sim_conf_filepath)

    # Return this object.
    return sim_state
