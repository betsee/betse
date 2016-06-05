#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Medium-level simulation configuration classes wrapping low-level YAML
dictionaries deserialized from disk.
'''

#FIXME: Refactor all functions defined by "betse.science.config.sim_config"
#into methods of the "SimConfigWrapper" class defined below.

# ....................{ IMPORTS                            }....................
# import yaml
# from betse import pathtree
from betse.science.config import sim_config
from betse.util.path import files, paths
# from betse.util.type import types
# from betse.util.io.log import logs

# ....................{ CLASSES                            }....................
class SimConfigWrapper(object):
    '''
    Medium-level simulation configuration wrapper wrapping a low-level
    dictionary deserialized from a YAML-formatted simulation configuration file.

    This wrapper superficially wraps this dictionary with convenience methods
    safely modifying this dictionary. This wrapper is lower level than the
    high-level `Parameters` simulation configuration, which transforms this
    dictionary into numerous high-level objects rather than merely wrapping this
    dictionary. While the `Parameters` configuration is principally used by
    backend simulation modelling in the `betse.science` package, this wrapper is
    principally used by frontend logic modifying simulation configurations on
    behalf of either interactive users (e.g., BETSE's GUI) or automated tests.

    Attributes
    ----------
    _config : dict
        Low-level dictionary deserialized from `_filename`.
    _basename : str
        Basename of `_filename`.
    _dirname : str
        Absolute or relative path of the directory containing `_filename`.
    _filename : str
        Absolute or relative path of the YAML-formatted simulation configuration
        file deserialized into `_config`.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, filename: str) -> None:
        '''
        Wrap the low-level dictionary deserialized from the passed
        YAML-formatted simulation configuration file.

        Parameters
        ----------
        filename : str
            Absolute or relative path of this file.
        '''

        # Classify the passed parameters.
        self._filename = filename

        # Absolute or relative path of the directory containing this file.
        self._dirname = paths.get_dirname(filename)

        # Basename of this file.
        self._basename = paths.get_basename(filename)

        # Deserialize this YAML file into a dictionary.
        self._config = sim_config.read(filename)


    @classmethod
    def wrap_new_default(cls, filename: str) -> None:
        '''
        Write the default YAML-formatted simulation configuration to the passed
        path, recursively copy all external resources (e.g., geometry masks)
        referenced and hence required by this configuration into this path's
        directory, and return an instance of this class encapsulating this
        configuration.

        This factory method creates a valid simulation configuration consumable
        by all BETSE CLI commands (e.g., `betse sim`), modified from the default
        simulation configuration shipped with BETSE as follows:

        * The `plot after solving` option in the `results options` section is
          coerced to `False`, preventing hapless end-users from drowning under
          an intimidating deluge of plot windows irrelevant to "beginner" usage.

        Parameters
        ----------
        filename : str
            Absolute or relative path of the simulation configuration file to be
            written. Since this file will be YAML-formatted, this filename
            should ideally be suffixed by a valid YAML filetype: namely, either
            `.yml` or `.yaml`. This is _not_ strictly necessary, but is strongly
            recommended.

        Raises
        ----------
        BetseExceptionFile
            If this file already exists.
        '''

        # Create this YAML file.
        sim_config.write_default(filename)

        # Create and return an instance of this class wrapping this file.
        return cls(filename=filename)

    # ..................{ PROPERTIES                         }..................
    # For safety, these properties lack setters and hence are read-only.

    @property
    def basename(self) -> str:
        '''
        Basename of the configuration file wrapped by this encapsulation object.
        '''

        return self._basename


    @property
    def dirname(self) -> str:
        '''
        Absolute or relative path of the directory containing the configuration
        file wrapped by this encapsulation object.
        '''

        return self._dirname


    @property
    def filename(self) -> str:
        '''
        Absolute or relative path of the configuration file wrapped by this
        encapsulation object.
        '''

        return self._filename

    # ..................{ GETTERS                            }..................

    # ..................{ SETTERS                            }..................

    # ..................{ WRITERS                            }..................
    def overwrite(self) -> None:
        '''
        Reserialize the current low-level configuration dictionary to the
        current configuration file.

        This method silently overwrites the contents of this file with the
        (possibly modified) contents of this dictionary.
        '''

        # Delete this configuration file, preventing the subsequent write from
        # raising an ignorable exception.
        files.remove_if_found(self._filename)

        # Recreate this configuration file.
        self.write(self._filename)


    def write(self, filename: str) -> None:
        '''
        Serialize the current low-level configuration dictionary to the passed
        simulation configuration file in YAML format.

        If this file already exists, an exception is raised.

        Parameters
        ----------
        filename : str
            Absolute or relative path of the simulation configuration file to be
            written. Since this file will be YAML-formatted, this filename
            should ideally be suffixed by a valid YAML filetype: namely, either
            `.yml` or `.yaml`. This is _not_ strictly necessary, but is strongly
            recommended.

        Raises
        ----------
        BetseExceptionFile
            If this file already exists.
        '''

        sim_config.write(filename, self._config)

    # ..................{ DISABLERS                          }..................
    #FIXME: The implementation of the following methods is fundamentally unsafe.
    #If the structure of the underlying YAML file changes, these methods could
    #silently fail (e.g., if the "plot while solving" option were renamed to
    #"is plotting during"). To combat this, all attempts to directly modify the
    #"self._config" dictionary below *MUST* instead defer to a newly defined
    #set_config_option() method accepting one or more key names followed by the
    #value to set such keys to: e.g.,
    #
    #    set_config_option(('results options', 'plot while solving'), False)
    #
    #If the passed configuration option does *NOT* already exist, that method
    #should raise a human-readable exception. Inevitable future problem solved!

    def disable_interaction(self) -> None:
        '''
        Disable all simulation configuration options either requiring
        interactive user input _or_ displaying graphical output intended for
        interactive user consumption (e.g., plots, animations).

        This method is intended to be called by non-interactive automation
        (e.g., tests, scripts) expecting simulations to behave silently.
        '''

        # Disable display of all mid- and post-simulation plots and animations.
        self._config['results options']['plot while solving'] = False
        self._config['results options']['plot after solving'] = False

    # ..................{ MINIMIZERS                         }..................
    def minify(self) -> None:
        '''
        Minimize the space and time costs associated with running the simulation
        configured by this configuration while preserving all fundamental
        configuration features.

        Specifically, this method numerically reduces all configuration options
        pertaining to either world size _or_ simulation time to their **minimum
        permissible values** (i.e., the smallest values still preserving
        simulation stability). This method is intended to be called only by test
        automation.
        '''

        # Minify initialization time to exactly ten sampled time steps. For
        # safety, permit currently defined options smaller than the minimums
        # defined below to override these minimums.
        init = self._config['init time settings']
        # Duration of each time step in seconds.
        init['time step'] = min(
            float(init['time step']), 1.0e-3)
        # Interval to sample such steps at in seconds.
        init['sampling rate'] = min(
            float(init['sampling rate']), init['time step'])
        # Total simulation time in seconds.
        init['total time'] = min(
            float(init['total time']), 10.0e-3)

        # Minify simulation time to the same durations.
        sim = self._config['sim time settings']
        sim['time step'] = min(
            float(sim['time step']), init['time step'])
        sim['sampling rate'] = min(
            float(sim['sampling rate']), sim['time step'])
        sim['total time'] = min(
            float(sim['total time']), init['total time'])

        # Minify the physical dimensions of the cell cluster in meters. By
        # experimentation, the default simulation configuration exhibits
        # instabilities (e.g., on performing the cutting event) raising fatal
        # exceptions for physical dimensions less than that specified below.
        # Hence, this appears to currently be a fairly hard minimum.
        world = self._config['world options']
        world['world size'] = min(float(world['world size']), 100e-6)

        # Minify ECM-specific grid size. For similar reasons as above, the
        # computational grid size specified below appears to be a hard minimum.
        ecm = self._config['general options']
        ecm['comp grid size'] = min(int(ecm['comp grid size']), 20)
        ecm['plot grid size'] = min(int(ecm['plot grid size']), 50)
