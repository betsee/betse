#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level simulation configuration classes wrapping low-level dictionaries both
serialized to and deserialized from on-disk YAML-formatted files.
'''

#FIXME: Refactor all functions defined by the "betse.science.config.confio"
#submodule into methods of the "SimConfigWrapper" class defined below.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseNumericException
from betse.science.config import confdefault, confio
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.type.types import type_check, NumericTypes

# ....................{ SUPERCLASSES                       }....................
class SimConfigWrapper(object):
    '''
    High-level simulation configuration wrapper wrapping a low-level
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
        self._config = confio.read(filename)


    #FIXME: Eliminate this method in favour of the existing
    #sim_config.write() function *AFTER* refactoring that method as
    #documented in that submodule.
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
        BetseFileException
            If this file already exists.
        '''

        # Create this YAML file.
        confdefault.write(filename)

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

    # ..................{ WRITERS                            }..................
    def overwrite(self) -> None:
        '''
        Reserialize the current low-level configuration dictionary to the
        current configuration file.

        This method silently overwrites the contents of this file with the
        (possibly modified) contents of this dictionary.
        '''

        # Log this overwrite attempt.
        logs.log_info('Overwriting current simulation configuration...')

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
        BetseFileException
            If this file already exists.
        '''

        confio.write(filename, self._config)

    # ..................{ PROPERTIES ~ bool                  }..................
    @property
    def is_ecm(self) -> bool:
        '''
        `True` only if the extracellular matrix (ECM) is simulated by this
        configuration.
        '''

        return self._config['general options']['simulate extracellular spaces']


    @is_ecm.setter
    @type_check
    def is_ecm(self, is_ecm: bool) -> None:
        '''
        Set whether the extracellular matrix (ECM) is simulated by this
        configuration.
        '''

        self._config['general options']['simulate extracellular spaces'] = (
            is_ecm)

    # ..................{ PROPERTIES ~ bool : networks       }..................
    @property
    def is_brn(self) -> bool:
        '''
        `True` only if the biochemical reaction network (BRN) is simulated by
        this configuration.
        '''

        return self._config['metabolism settings']['metabolism simulated']


    @property
    def is_grn(self) -> bool:
        '''
        `True` only if the gene regulatory network (GRN) is simulated by
        this configuration.
        '''

        return self._config['gene regulatory network settings']\
            ['gene regulatory network simulated']


    @is_brn.setter
    @type_check
    def is_brn(self, is_brn: bool) -> None:
        '''
        Set whether the biochemical reaction network (BRN) is simulated by this
        configuration.
        '''

        self._config['metabolism settings']['metabolism simulated'] = is_brn


    @is_grn.setter
    @type_check
    def is_grn(self, is_grn: bool) -> None:
        '''
        Set whether the gene regulatory network (GRN) is simulated by this
        configuration.
        '''

        self._config['gene regulatory network settings']\
            ['gene regulatory network simulated'] = is_grn

    # ..................{ PROPERTIES ~ bool : plot           }..................
    @property
    def is_plot_after_sim(self) -> bool:
        '''
        `True` only if post-simulation plots are enabled by this configuration.
        '''

        return self._config['results options']\
            ['after solving']['plots']['enabled']


    @is_plot_after_sim.setter
    @type_check
    def is_plot_after_sim(self, is_plot_after_sim: bool) -> None:
        '''
        Set whether post-simulation plots are enabled by this configuration.
        '''

        self._config['results options']\
            ['after solving']['plots']['enabled'] = is_plot_after_sim

    # ..................{ PROPERTIES ~ bool : anim : getter  }..................
    @property
    def is_anim_while_sim(self) -> bool:
        '''
        `True` only if in-simulation animations are enabled by this
        configuration.
        '''

        return self._config['results options']\
            ['while solving']['animations']['enabled']


    @property
    def is_anim_while_sim_save(self) -> bool:
        '''
        `True` only if in-simulation animations are saved by this configuration.
        '''

        return self._config['results options']\
            ['while solving']['animations']['save']


    @property
    def is_anim_after_sim(self) -> bool:
        '''
        `True` only if post-simulation animations are enabled by this
        configuration.
        '''

        return self._config['results options']\
            ['after solving']['animations']['enabled']


    @property
    def is_anim_after_sim_save(self) -> bool:
        '''
        `True` only if post-simulation animations are saved by this
        configuration.
        '''

        return self._config['results options']\
            ['after solving']['animations']['save']

    # ..................{ PROPERTIES ~ bool : anim : setter  }..................
    @is_anim_while_sim.setter
    @type_check
    def is_anim_while_sim(self, is_anim_while_sim: bool) -> None:
        '''
        Set whether in-simulation animations are enabled by this configuration.
        '''

        self._config['results options']\
            ['while solving']['animations']['enabled'] = is_anim_while_sim


    @is_anim_while_sim_save.setter
    @type_check
    def is_anim_while_sim_save(self, is_anim_while_sim_save: bool) -> None:
        '''
        Set whether in-simulation animations are saved by this configuration.
        '''

        self._config['results options']\
            ['while solving']['animations']['save'] = is_anim_while_sim_save


    @is_anim_after_sim.setter
    @type_check
    def is_anim_after_sim(self, is_anim_after_sim: bool) -> None:
        '''
        Set whether post-simulation animations are enabled by this
        configuration.
        '''

        self._config['results options']\
            ['after solving']['animations']['enabled'] = is_anim_after_sim


    @is_anim_after_sim_save.setter
    @type_check
    def is_anim_after_sim_save(self, is_anim_after_sim_save: bool) -> None:
        '''
        Set whether post-simulation animations are saved by this configuration.
        '''

        self._config['results options']\
            ['after solving']['animations']['save'] = is_anim_after_sim_save

    # ..................{ PROPERTIES ~ float                 }..................
    @property
    def environment_size(self) -> float:
        '''
        Square dimension in meters of the external environment containing the
        cell cluster for this configuration.

        For simplicity, BETSE constrains the environment to be square in shape.
        This dimension thus defines both the environmental width _and_ height.
        '''

        # Coerce the current number to a float for safety.
        return float(self._config['world options']['world size'])


    @environment_size.setter
    @type_check
    def environment_size(self, environment_size: float) -> None:
        '''
        Set the square dimension in meters of the external environment
        containing the cell cluster for this configuration.
        '''

        # Coerce the passed number to a float for safety.
        self._config['world options']['world size'] = float(environment_size)

    # ..................{ PROPERTIES ~ time : sim : getter   }..................
    @property
    def sim_sample_rate(self) -> float:
        '''
        Simulation-specific sample rate in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._config['sim time settings']['sampling rate'])


    @property
    def sim_time(self) -> float:
        '''
        Simulation-specific duration in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._config['sim time settings']['total time'])


    @property
    def sim_time_step(self) -> float:
        '''
        Simulation-specific time step in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._config['sim time settings']['time step'])

    # ..................{ PROPERTIES ~ time : sim : setter   }..................
    @sim_sample_rate.setter
    @type_check
    def sim_sample_rate(self, sim_sample_rate: NumericTypes) -> None:
        '''
        Set the simulation-specific time step in seconds for this configuration
        to the passed integer or float.
        '''

        # If this sample rate is less than the current time step, raise an
        # exception.
        if sim_sample_rate < self.sim_time_step:
            raise BetseNumericException(
                'Simulation sample rate {} < time step {}.'.format(
                    sim_sample_rate, self.sim_time_step))

        # Coerce the passed number to a float for safety.
        self._config['sim time settings']['sampling rate'] = float(
            sim_sample_rate)


    @sim_time.setter
    @type_check
    def sim_time(self, sim_time: NumericTypes) -> None:
        '''
        Set the simulation-specific duration in seconds for this configuration
        to the passed integer or float.
        '''

        # If this duration is less than the current sample rate, raise an
        # exception.
        if sim_time < self.sim_sample_rate:
            raise BetseNumericException(
                'Simulation duration {} < sample rate {}.'.format(
                    sim_time, self.sim_sample_rate))

        # Coerce the passed number to a float for safety.
        self._config['sim time settings']['total time'] = float(sim_time)


    @sim_time_step.setter
    @type_check
    def sim_time_step(self, sim_time_step: NumericTypes) -> None:
        '''
        Set the simulation-specific time step in seconds for this configuration
        to the passed integer or float.
        '''

        # If this time step is non-positive, raise an exception.
        if sim_time_step <= 0:
            raise BetseNumericException(
                'Simulation time step {} <= 0.'.format(sim_time_step))

        # Coerce the passed number to a float for safety.
        self._config['sim time settings']['time step'] = float(sim_time_step)

    # ..................{ PROPERTIES ~ ion profile           }..................
    @property
    def ion_profile(self) -> str:
        '''
        Name of the currently enabled ion profile for this configuration.
        '''

        return self._config['general options']['ion profile']


    #FIXME: Validate the passed string. Presumably, we have logic elsewhere
    #already doing so. Leverage such logic here.
    @ion_profile.setter
    @type_check
    def ion_profile(self, ion_profile_name: str) -> None:
        '''
        Set the ion profile for this configuration to that with the passed name.

        Parameters
        ----------
        ion_profile_name : str
            Name of the ion profile to be enabled. See the default simulation
            configuration file for a list of all supported strings.
        '''

        self._config['general options']['ion profile'] = ion_profile_name
