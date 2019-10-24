#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixture classes encapsulating test-related simulation configurations.
'''

#FIXME: Most use of the increasingly obsolete "SimConfTestInternal.config"
#wrapper attribute (both here and everywhere else) should be replaced by use of
#the new "SimConfTestInternal.p" property, which increasingly provides all test
#functionality.

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.types import type_check
from py._path.local import LocalPath

# ....................{ SUPERCLASSES                      }....................
class SimConfTestABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all simulation configuration context subclasses,
    each encapsulating state and metadata for a test-related simulation
    configuration.

    Simulation configuration fixtures typically return instances of this class
    as a means of communicating this context to other fixtures and tests.

    Attributes
    ----------
    conf_dirname : str
        Absolute dirname of the directory of the simulation configuration file
        wrapped by this wrapper (i.e., of :attr:`conf_filename`).
    conf_filename : str
        Absolute filename of the simulation configuration file wrapped by this
        wrapper.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official :mod:`py.path` submodule documentation.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self, conf_filename: str) -> None:
        '''
        Initialize this simulation configuration context.

        Parameters
        ----------
        conf_filename : str
            Absolute filename of the simulation configuration file encapsulated
            by this context. Although this file need *not* physically exist
            before this object is instantiated (i.e., before this method is
            called), the subclass should ensure this file exists once this
            object is instantiated (i.e., immediately before the subclass
            implementation of this method returns).
        '''

        # Defer heavyweight imports.
        from betse.util.path import pathnames

        # Classify all passed parameters publicly, permitting access by tests.
        self.conf_filename = conf_filename
        self.conf_dirname = pathnames.get_dirname(self.conf_filename)

    # ..................{ CONTEXTS                          }..................
    def context(self) -> 'contextlib.contextmanager':
        '''
        Context manager changing the current working directory (CWD) of the
        current test to the directory containing this configuration file for
        the duration of this context.

        Default simulation configuration paths are relative to the directory
        containing the simulation configuration file: namely, this temporary
        directory. Changing directories resolves these paths to this directory.
        (Failing to do so would incorrectly resolve these paths to the current
        directory with predictably disastrous outcomes.) While this class could
        instead globally search-and-replace all relative simulation
        configuration paths with absolute paths, doing so would be considerably
        more complex, fragile, and error-prone than changing directories.
        '''

        # Defer heavyweight imports.
        from betse.util.os.shell import shelldir

        # Defer to the generator returned by the following utility function.
        return shelldir.setting_cwd(self.conf_dirname)

    # ..................{ SUBCLASS                          }..................
    # Subclasses are required to define the following read-only properties.

    @abstractproperty
    def p(self) -> 'betse.science.parameters.Parameters':
        '''
        High-level simulation configuration encapsulated by this test wrapper.
        '''

        pass

# ....................{ SUBCLASSES                        }....................
class SimConfTestInternal(SimConfTestABC):
    '''
    Simulation configuration context subclass encapsulating a temporary
    simulation configuration file internally isolated to the current test.

    Since this configuration is guaranteed to be test-specific and hence safely
    modifiable, caller fixtures and tests are advised to modify the contents of
    this configuration (e.g., to exercise specific feature sets and edge
    cases).

    Attributes
    ----------
    config : SimConfigWrapper
        Simulation configuration wrapper wrapping the low-level dictionary
        deserialized from the YAML-formatted simulation configuration file with
        path :attr:`conf_filepath`. Note the contents of this dictionary may be
        desynchronized from those of this file. For efficiency, callers may
        modify this dictionary to suite test requirements *before*
        reserializing this dictionary back to this file.
    conf_filepath : LocalPath
        Absolute filename of a temporary simulation configuration file specific
        to the parent fixture as a :class:`py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented :mod:`os.path`
        module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official :mod:`py.path` submodule documentation.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,
        src_conf_filename: str,
        trg_conf_filepath: LocalPath,
    ) -> None:
        '''
        Initialize this context by copying the passed source to target
        simulation configuration file.

        Specifically, this method (in order):

        #. Copies the passed source to target simulation configuration file.
        #. Copies *all* external resources (e.g., image files) referenced and
           hence required by the passed source simulation configuration file to
           the directory of the passed target simulation configuration file.
        #. Sanitizes the copied simulation configuration file for all child
           fixtures and tests by unconditionally disabling options either
           requiring interactive input *or* displaying interactive output.

        This method does *not* minimize the space and time costs of running
        simulations configured by this configuration. Doing so is typically
        desirable but can obscure edge-case issues, including computational
        instability produced by the default non-minified time steps.

        Caveats
        ----------
        For efficiency, note that the configuration changes applied by this
        method (listed above) reside *only* in-memory; they have yet to be
        written back to disk. Callers are required to do so manually (e.g., via
        the optional ``is_overwrite_conf`` parameter passed to the
        :meth:`CLISimTester.run_subcommands` method).

        Parameters
        ----------
        src_conf_filename : str
            Absolute filename of the source simulation configuration file from
            which this method safely copies the passed target simulation
            configuration file to.
        trg_conf_filepath : LocalPath
            Absolute filename of the target simulation configuration file to
            which this method safely copies both the passed source simulation
            configuration file *and* all external resources (e.g., image files)
            referenced and hence required by that file. For compliance with
            :mod:`pytest` internals, this filename is a high-level
            :class:`py.path.local` instance rather than a low-level string.

        Raises
        ----------
        BetseDirException
            If any external file required by this target file already exists.
        BetseFileException
            If this target file already exists.
        '''

        # Defer heavyweight imports. Notably, the "simconfwrapper" submodule
        # inherits from the main codebase and is thus *NOT* safely importable
        # at module scope.
        from betse.science.parameters import Parameters
        from betse_test.fixture.simconf.simconfwrapper import (
            SimConfigTestWrapper)

        # Absolute filename of this file.
        trg_conf_filename = str(trg_conf_filepath)

        # Initialize our superclass with this filename.
        super().__init__(conf_filename=trg_conf_filename)

        # Classify the passed parameters. While the "self.config" object
        # classified below provides this filename as a low-level string, this
        # high-level "py.path.local" instance is useful in fixtures and tests.
        self.conf_filepath = trg_conf_filepath

        # Copy this source to target file.
        p = Parameters()
        p.copy(
            src_conf_filename=src_conf_filename,
            trg_conf_filename=trg_conf_filename)

        # Test-specific wrapper encapsulating this file.
        self.config = SimConfigTestWrapper(p)

        # Sanitize this configuration for all child fixtures and tests.
        self.config.disable_interaction()

    # ..................{ SUPERCLASS                        }..................
    @property
    def p(self) -> 'betse.science.parameters.Parameters':
        return self.config.p


class SimConfTestExternal(SimConfTestABC):
    '''
    Simulation configuration context subclass encapsulating an external
    simulation configuration file residing at any directory and hence *not*
    necessarily isolated to the current test.

    Attributes
    ----------
    _p : Parameters
        High-level simulation configuration encapsulating this external
        simulation configuration file.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize this simulation configuration context.
        '''

        # Defer heavyweight imports.
        from betse.science.parameters import Parameters

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # In-memory simulation configuration deserialized from this file.
        self._p = Parameters.make(conf_filename=self.conf_filename)

    # ..................{ SUPERCLASS                        }..................
    @property
    def p(self) -> 'betse.science.parameters.Parameters':
        return self._p
