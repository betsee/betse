#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: This submodule requires heavy refactoring away from the current
#low-level approach in favour of the new high-level "confabc"-based approach --
#namely, the "SimConfABC" base class coupled with the conf_alias() data
#descriptor. Together, this new functionality provides a substantially better
#YAML-to-Python-object-mapping (YPOM) than the ad-hoc and overly verbose
#boilerplate implemented below. Specifically:
#
#* Eliminate all properties defined below.
#* Replace all usage of the low-level "self._p._conf" dictionary with high-level
#  data descriptors defined by "self._p".

'''
Test-specific simulation configuration classes wrapping low-level dictionaries
both serialized to and deserialized from on-disk YAML-formatted files.
'''

# ....................{ IMPORTS                            }....................
# This subclass necessarily imports from submodules defined by the main codebase
# and is thus *NOT* safely importable in a fixture submodule directly imported
# by a "conftest" plugin module. To defer the importation of this submodule
# until *AFTER* test collection, this submodule is intentionally segregated.

from betse.exceptions import BetseNumericException
from betse.science.config import confdefault, confio
from betse.science.parameters import Parameters
from betse.science.simulate.pipe import piperunreq
from betse.science.visual.anim.animpipe import AnimCellsPipe
from betse.science.visual.plot.pipe.plotpipecell import PlotCellPipe
from betse.science.visual.plot.pipe.plotpipecells import PlotCellsPipe
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.type.types import type_check, NumericTypes

# ....................{ SUPERCLASSES                       }....................
class SimConfigTestWrapper(object):
    '''
    Test-specific simulation configuration wrapper wrapping a low-level
    dictionary deserialized from a YAML-formatted simulation configuration file.

    This wrapper is intended to be instantiated *only* by non-interactive
    automation (e.g., tests, scripts).

    This wrapper superficially wraps this dictionary with convenience methods
    safely modifying this dictionary. This wrapper is lower level than the
    high-level :class:`Parameters` simulation configuration, which transforms
    this dictionary into numerous high-level objects rather than merely wrapping
    this dictionary. While the :class:`Parameters` configuration is principally
    used by backend simulation modelling in the `betse.science` package, this
    wrapper is principally used by frontend logic modifying simulation
    configurations on behalf of either interactive users (e.g., BETSE's GUI) or
    automated tests.

    Attributes
    ----------
    _config : dict
        Low-level dictionary deserialized from :attr:`_filename`.
    _basename : str
        Basename of :attr:`_filename`.
    _dirname : str
        Absolute or relative path of the directory containing :attr:`_filename`.
    _filename : str
        Absolute or relative path of the YAML-formatted simulation configuration
        file deserialized into :attr:`_config`.
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

        # Deserialize this file into a high-level in-memory object.
        self._p = Parameters(config_filename=filename)


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


    @type_check
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

        confio.write(filename, self._p._conf)

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
        return float(self._p._conf['world options']['world size'])


    @environment_size.setter
    @type_check
    def environment_size(self, environment_size: float) -> None:
        '''
        Set the square dimension in meters of the external environment
        containing the cell cluster for this configuration.
        '''

        # Coerce the passed number to a float for safety.
        self._p._conf['world options']['world size'] = float(environment_size)

    # ..................{ PROPERTIES ~ time : sim : getter   }..................
    @property
    def sim_sample_rate(self) -> float:
        '''
        Simulation-specific sample rate in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._p._conf['sim time settings']['sampling rate'])


    @property
    def sim_time(self) -> float:
        '''
        Simulation-specific duration in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._p._conf['sim time settings']['total time'])


    @property
    def sim_time_step(self) -> float:
        '''
        Simulation-specific time step in seconds for this configuration.
        '''

        # Coerce the current number to a float for safety.
        return float(self._p._conf['sim time settings']['time step'])

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
        self._p._conf['sim time settings']['sampling rate'] = float(
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
        self._p._conf['sim time settings']['total time'] = float(sim_time)


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
        self._p._conf['sim time settings']['time step'] = float(sim_time_step)

    # ..................{ PROPERTIES ~ ion profile           }..................
    @property
    def ion_profile(self) -> str:
        '''
        Name of the currently enabled ion profile for this configuration.
        '''

        return self._p._conf['general options']['ion profile']


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

        self._p._conf['general options']['ion profile'] = ion_profile_name

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
        init = self._p._conf['init time settings']
        # Duration of each time step in seconds.
        init['time step'] = min(float(init['time step']), 1.0e-3)
        # Interval to sample such steps at in seconds.
        init['sampling rate'] = min(
            float(init['sampling rate']), init['time step'])
        # Total simulation time in seconds. The first digit effectively defines
        # the number of sampled time steps, by the above choice of time step.
        init['total time'] = min(float(init['total time']), 3.0e-3)

        # Minify simulation time to the same durations. To ensure that
        # property setter validation compares durations in the expected manner,
        # properties are assigned in order of increasing duration.
        self.sim_time_step =   min(self.sim_time_step, init['time step'])
        self.sim_sample_rate = min(self.sim_sample_rate, self.sim_time_step)
        self.sim_time =        min(self.sim_time, init['total time'])

        # Minify the physical dimensions of the cell cluster in meters. By
        # experimentation, the default simulation configuration both:
        #
        # * Exhibits expected instabilities for physical dimensions less than
        #   100e-6, typically due to the cell cluster cutting event. Since these
        #   instabilities are expected and hence best ignored, the dimensions
        #   specified below should *NEVER* be smaller than 100e-6.
        # * Exhibits unexpected instabilities for physical dimensions greater
        #   than 250e-6, typically due to non-linear interactions between
        #   channel electroosmosis and physical deformation. Since these
        #   instabilities are unexpected and hence *NOT* safely ignorable, the
        #   dimensions specified below should *ALWAYS* be in this range.
        #
        # Hence, these dimensions reside in an optimal range both avoiding
        # expected instabilities *AND* exposing unexpected instabilities.
        self.environment_size = min(self.environment_size, 100e-6)

        #FIXME: Uncommenting the following line reliably causes the following
        #test to fail with a computational instability during the "sim" phase:
        #
        #    ./test -k test_cli_sim_ecm
        #
        #Since it's unclear whether this behaviour is expected or not, this line
        #remains commented out. (If this behaviour is indeed unexpected, this
        #line should be shifted into the enable_exports_all() method and the
        #more preferable global default retained above).

        # self.environment_size = min(self.environment_size, 250e-6)

        # Minify ECM-specific grid size. For similar reasons as above, the
        # computational grid size specified below appears to be a hard minimum.
        ecm = self._p._conf['general options']
        ecm['comp grid size'] = min(int(ecm['comp grid size']), 20)
        ecm['plot grid size'] = min(int(ecm['plot grid size']), 50)

    # ..................{ DISABLERS                          }..................
    #FIXME: The implementation of the following methods is fundamentally unsafe.
    #If the structure of the underlying YAML file changes, these methods could
    #silently fail (e.g., if the "plot while solving" option were renamed to
    #"is plotting during"). To combat this, all attempts to directly modify the
    #"self._p._conf" dictionary below *MUST* instead defer to a newly defined
    #set_config_option() method accepting one or more key names followed by the
    #value to set such keys to: e.g.,
    #
    #    set_config_option(('results options', 'plot while solving'), False)
    #
    #If the passed configuration option does *NOT* already exist, that method
    #should raise a human-readable exception. Inevitable future problem solved!

    def disable_visuals(self) -> None:
        '''
        Disable all visual exports, including displaying and saving of all in-
        and post-simulation plots and animations.
        '''

        self._p.anim.is_after_sim = False
        self._p.anim.is_while_sim = False
        self._p.plot.is_after_sim = False


    def disable_interaction(self) -> None:
        '''
        Disable all simulation configuration options either requiring
        interactive user input _or_ displaying graphical output intended for
        interactive user consumption (e.g., plots, animations).

        This method is intended to be called by non-interactive automation
        (e.g., tests, scripts) expecting simulations to behave silently.
        '''

        results = self._p._conf['results options']
        results['while solving']['animations']['show'] = False
        results['after solving']['plots']['show'] = False
        results['after solving']['animations']['show'] = False

    # ..................{ ENABLERS                           }..................
    def enable_networks(self) -> None:
        '''
        Enable both biochemical reaction and gene regulatory networks.
        '''

        self._p._conf['metabolism settings']['metabolism simulated'] = True
        self._p._conf['gene regulatory network settings'][
            'gene regulatory network simulated'] = True


    def enable_vg_ion_channels_all(self) -> None:
        '''
        Enable all voltage-gated ion channels (e.g., sodium, potassium) _and_
        all features required by these channels.

        This method is intended to be called by non-interactive test suites
        exercising these channels. Specifically, this method enables:

        * The extracellular matrix (ECM).
        * The mammalian ion profile (i.e., `animal`), enabling all ions.
        * The intervention increasing the permeability of all cell membranes to
          sodium (Na+).
        * The voltage-gated sodium (Na+) channel `Nav1p2`, corresponding to the
          adult human brain.
        * The voltage-gated potassium (K+) channel `K_Slow`.
        * Decreased time step and sampling rates, ensuring simulation stability.
        * Increased duration and cell count, exposing simulation instabilities.

        For efficiency, this method disables all visuals -- including both in-
        and post-simulation animations and plots.
        '''

        # For efficiency, disable all visuals.
        self.disable_visuals()

        # Enable all features required by these channels.
        self._p.sim_ECM = True
        self.ion_profile = 'animal'

        # For stability, decrease both the time step and sampling rates.
        self.sim_time_step =   1e-4
        self.sim_sample_rate = 1e-3

        # For completeness, increase both the duration and cell count.
        self.sim_time = 50e-3
        self.environment_size = 250e-6

        # Enable the intervention increasing sodium membrane permeability.
        # Although the current default values for this intervention track those
        # defined below fairly closely, the latter are nonetheless explicitly
        # defined below to avoid issues when the former inevitably change.
        sodium_membrane_permeability = self._p._conf['change Na mem']
        sodium_membrane_permeability['event happens'] = True
        sodium_membrane_permeability['change rate']   =  1.0e-3
        sodium_membrane_permeability['change start']  =  5.0e-3
        sodium_membrane_permeability['change finish'] = 30.0e-3
        sodium_membrane_permeability['apply to'] = ['spot',]

        # Enable the voltage-gated sodium (Na+) channel Nav1p2.
        voltage_gated_sodium_channel = self._p._conf['voltage gated Na+']
        voltage_gated_sodium_channel['turn on'] = True
        voltage_gated_sodium_channel['channel type'] = ['Nav1p2',]
        # voltage_gated_sodium_channel['max value'] = 5.0e-6
        voltage_gated_sodium_channel['apply to'] = ['all',]

        # Enable the voltage-gated potassium (K+) channel K_Slow.
        voltage_gated_potassium_channel = self._p._conf['voltage gated K+']
        voltage_gated_potassium_channel['turn on'] = True
        voltage_gated_potassium_channel['channel type'] = ['K_Slow',]
        # voltage_gated_potassium_channel['max value'] = 5.0e-7
        voltage_gated_potassium_channel['apply to'] = ['all',]

    # ..................{ ENABLERS ~ export                  }..................
    @type_check
    def enable_anim_video(self, writer_name: str, filetype: str) -> None:
        '''
        Enable encoding of all enabled animations as compressed video of the
        passed filetype with the preferred matplotlib animation writer of the
        passed name.

        Parameters
        ----------
        writer_name : str
            Name of the matplotlib animation writer with which to encode video
            (e.g., `ffmpeg`, `imagemagick`).
        filetype : str
            Filetype of videos to encode with this writer (e.g., `mkv`, `mp4`).
        '''

        # Enable animations and animation saving in the general sense.
        self.enable_visuals_save()

        # Localize nested dictionaries for convenience.
        video = self._p._conf['results options']['save']['animations']['video']

        # Enable encoding of the passed filetype with the passed writer type.
        # For determinism, mandate that *ONLY* this writer (rather than two or
        # more writers) be used to do so.
        video['enabled'] = True
        video['filetype'] = filetype
        video['writers'] = [writer_name,]


    def enable_visuals_save(self) -> None:
        '''
        Enable saving of all visual exports, including all in- and
        post-simulation plots and animations.

        This method does *not* enable specific visual exports or simulation
        features required by specific visual exports.
        '''

        self._p.anim.is_while_sim_save = True
        self._p.anim.is_after_sim_save = True
        self._p.plot.is_after_sim_save = True

    # ..................{ ENABLERS ~ export : ecm            }..................
    def enable_exports_ecmless(self) -> None:
        '''
        Enable all available exports (e.g., CSVs, plots, animations) excluding
        those requiring extracellular spaces, all features required by these
        plots and animations, and any additional features trivially enabled
        *without* substantially increasing time or space complexity.

        Specifically, this method enables all exports, features, and settings
        enabled by the :meth:`_enable_visuals_common` method.
        '''

        # Enable all optional general settings supported by visual exports.
        self._enable_visuals_common()

        # Disable extracellular spaces.
        self._p.sim_ECM = False

        # For each type of export pipeline to be exercised and list of all
        # currently enabled exporters in this pipeline...
        for pipe_type, pipe_list in self._pipes_type_list:
            # For the name and metadata of each exporter (enabled or not)
            # supported by this pipeline...
            for pipe_exporter_name, pipe_exporter in pipe_type.iter_runners():
                # If this export needs extracellular spaces, ignore this export.
                if piperunreq.ECM in pipe_exporter.requirements:
                    continue
                # Else, this export does *NOT* need extracellular space.

                # New default export of this type appended to this pipeline.
                pipe_exporter_conf = pipe_list.append_default()
                pipe_exporter_conf.name = pipe_exporter_name


    def enable_exports_ecm(self) -> None:
        '''
        Enable all available exports (e.g., CSVs, plots, animations) including
        those requiring extracellular spaces, all features required by these
        plots and animations, and any additional features trivially enabled
        *without* substantially increasing time or space complexity.

        Specifically, this method enables:

        * All exports, features, and settings enabled by the
          :meth:`_enable_visuals_common` method.
        * The ``current_total``, ``electric_total``, ``fluid_total``, and
          ``voltage_total`` plots and animations by enabling:
          * The extracellular matrix (ECM).
        '''

        # Enable all optional general settings supported by visual exports.
        self._enable_visuals_common()

        # Enable extracellular spaces.
        self._p.sim_ECM = True

        # For each type of export pipeline to be exercised and list of all
        # currently enabled exporters in this pipeline...
        for pipe_type, pipe_list in self._pipes_type_list:
            # For the name of each exporter supported by this pipeline...
            for pipe_exporter_name, _ in pipe_type.iter_runners():
                # New default export of this type appended to this pipeline.
                pipe_exporter_conf = pipe_list.append_default()
                pipe_exporter_conf.name = pipe_exporter_name


    def _enable_visuals_common(self) -> None:
        '''
        Enable all visual exports (e.g., in- and post-simulation plots and
        animations) and simulation features required by these exports *without*
        enabling extracellular spaces.

        This method additionally enables optional settings improving test
        coverage but *not* explicitly required by these exports. Specifically,
        this method enables:

        * Saving of all visual exports.
        * Cell enumeration, labelling each cell by its 0-based index.
        * Current overlays, displaying current density streamlines.
        * All ion concentration plots and animations (e.g., calcium (Ca+),
          hydrogen (H+; pH)) by enabling:
          * The mammalian ion profile (i.e., ``animal``), enabling all ions.
        * The deformation plot and animation by enabling:
          * Galvanotaxis (i.e., deformations).
        * The ``fluid_intra`` plot and animation by enabling:
          * Fluid flow.
        * The ``pressure_mechanical`` plot and animation by enabling:
          * The mechanical pressure event.
        * The ``pressure_osmotic`` plot and animation by enabling:
          * Osmotic pressure.
        * The ``pump_density`` plot and animation by enabling:
          * Channel electroosmosis.
        * The ``voltage_polarity`` plot and animation by enabling:
          * Cell polarizability.
        * All other plots and animations *not* requiring extracellular spaces.
        '''

        # Enable saving of these exports.
        self.enable_visuals_save()

        # Localize nested dictionaries for convenience.
        results = self._p._conf['results options']
        variable = self._p._conf['variable settings']

        # Enable all simulation features required by these exports.
        self.ion_profile = 'animal'
        self._p._conf['apply pressure']['event happens'] = True
        variable['channel electroosmosis']['turn on'] = True
        variable['deformation']['turn on'] = True
        variable['fluid flow']['include fluid flow'] = True
        variable['pressures']['include osmotic pressure'] = True
        self._p.cell_polarizability = 1e-4

        # Enable all optional settings supported by these exports.
        results['enumerate cells'] = True
        results['overlay currents'] = True

    # ..................{ PRIVOTE ~ iterators                }..................
    @property
    def _pipes_type_list(self) -> tuple:
        '''
        Tuple of 2-tuples ``(pipe_type, pipe_list)``, describing each export
        pipeline to be exercised, where:

        * ``pipe_type`` is an instance of :class:`SimPipeABC`.
        * ``pipe_list`` is an instance of :class:`SimConfList` listing all
          currently enabled exporters in this pipeline.
        '''

        return (
            (PlotCellPipe,  self._p.plot.after_sim_pipeline_cell),
            (PlotCellsPipe, self._p.plot.after_sim_pipeline_cells),
            (AnimCellsPipe, self._p.anim.after_sim_pipeline),
        )
