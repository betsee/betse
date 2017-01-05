#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Test-specific simulation configuration classes wrapping low-level dictionaries
both serialized to and deserialized from on-disk YAML-formatted files.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.confwrap import SimConfigWrapper
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
# This subclass inherits a class defined by the main codebase and is thus *NOT*
# safely importable in a fixture submodule directly imported by a "conftest"
# plugin module. To permit the importation of this class to be deferred until
# *AFTER* test collection, this class is segregated into a separate submodule.

class SimConfigTestWrapper(SimConfigWrapper):
    '''
    Test-specific simulation configuration wrapper wrapping a low-level
    dictionary deserialized from a YAML-formatted simulation configuration file.

    This wrapper is intended to be instantiated _only_ by non-interactive
    automation (e.g., tests, scripts).
    '''

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
        # experimentation, the default simulation configuration exhibits
        # instabilities (e.g., on performing the cutting event) raising fatal
        # exceptions for physical dimensions less than that specified below.
        # Hence, this appears to currently be a fairly hard minimum.
        self.environment_size = min(self.environment_size, 100e-6)

        # Minify ECM-specific grid size. For similar reasons as above, the
        # computational grid size specified below appears to be a hard minimum.
        ecm = self._config['general options']
        ecm['comp grid size'] = min(int(ecm['comp grid size']), 20)
        ecm['plot grid size'] = min(int(ecm['plot grid size']), 50)

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

    def disable_visuals(self) -> None:
        '''
        Disable all visuals.

        Specifically, this method disables both display and export of all:

        * In-simulation animations.
        * Post-simulation animations.
        * Post-simulation plots.
        '''

        self.is_anim_while_sim = False
        self.is_anim_after_sim = False
        self.is_plot_after_sim = False


    def disable_interaction(self) -> None:
        '''
        Disable all simulation configuration options either requiring
        interactive user input _or_ displaying graphical output intended for
        interactive user consumption (e.g., plots, animations).

        This method is intended to be called by non-interactive automation
        (e.g., tests, scripts) expecting simulations to behave silently.
        '''

        results = self._config['results options']
        results['while solving']['animations']['show'] = False
        results['after solving']['plots']['show'] = False
        results['after solving']['animations']['show'] = False

    # ..................{ ENABLERS                           }..................
    def enable_anim_save(self) -> None:
        '''
        Enable both mid- and post-simulation animations _and_ the saving of
        these animations to disk in a general-purpose manner.

        This method does _not_ enable specific animations or features required
        by specific animations.
        '''

        self.is_anim_while_sim = True
        self.is_anim_after_sim = True
        self.is_anim_while_sim_save = True
        self.is_anim_after_sim_save = True


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
        self.enable_anim_save()

        # Localize nested dictionaries for convenience.
        video = self._config['results options']['save']['animations']['video']

        # Enable encoding of the passed filetype with the passed writer type.
        # For determinism, mandate that *ONLY* this writer (rather than two or
        # more writers) be used to do so.
        video['enabled'] = True
        video['filetype'] = filetype
        video['writers'] = [writer_name,]


    def enable_networks(self) -> None:
        '''
        Enable both biochemical reaction and gene regulatory networks.
        '''

        self.is_brn = True
        self.is_grn = True


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
        self.is_ecm = True
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
        sodium_membrane_permeability = self._config['change Na mem']
        sodium_membrane_permeability['event happens'] = True
        sodium_membrane_permeability['change rate']   =  1.0e-3
        sodium_membrane_permeability['change start']  =  5.0e-3
        sodium_membrane_permeability['change finish'] = 30.0e-3
        sodium_membrane_permeability['apply to'] = ['spot',]

        # Enable the voltage-gated sodium (Na+) channel Nav1p2.
        voltage_gated_sodium_channel = self._config['voltage gated Na+']
        voltage_gated_sodium_channel['turn on'] = True
        voltage_gated_sodium_channel['channel type'] = ['Nav1p2',]
        # voltage_gated_sodium_channel['max value'] = 5.0e-6
        voltage_gated_sodium_channel['apply to'] = ['all',]

        # Enable the voltage-gated potassium (K+) channel K_Slow.
        voltage_gated_potassium_channel = self._config['voltage gated K+']
        voltage_gated_potassium_channel['turn on'] = True
        voltage_gated_potassium_channel['channel type'] = ['K_Slow',]
        # voltage_gated_potassium_channel['max value'] = 5.0e-7
        voltage_gated_potassium_channel['apply to'] = ['all',]


    def enable_visuals_all(self) -> None:
        '''
        Enable all supported plots and animations, all features required by
        these plots and animations, and any additional features trivially
        enabled _without_ substantially increasing time or space complexity.

        This method is intended to be called by non-interactive test suites
        exercising all plots and animations. Specifically, this method enables:

        * Cell enumeration, labelling each cell by its 0-based index.
        * Current overlays, displaying current density streamlines.
        * The calcium (Ca) plot and animation by enabling:
          * The mammalian ion profile (i.e., `animal`), enabling all ions
            including calcium.
        * The deformation plot and animation by enabling:
          * Galvanotaxis (i.e., deformations).
        * The pH plot and animation dependent upon hydrogen (H) by enabling:
          * The mammalian ion profile (i.e., `animal`), enabling all ions
            including hydrogen.
        * The "Membrane" plot and animation of membrane pump density by
          enabling:
          * Membrane pump/channel movement via electrophoresis/osmosis.
        * The "P cell", "Osmotic P", and "Force" plots and animations of
          electroosmotic pressure, osmotic pressure, and hydrostatic body force
          respectively by enabling:
          * Osmotic pressure.
        * The "Vcell", "Venv", and "Current" plots and animations of cellular
          voltage, environmental voltage, and extracellular current by enabling:
          * The extracellular matrix (ECM).
        * The "Velocity" plot and animation of cluster fluid velocity by
          enabling:
          * Fluid flow.
          * Electrostatic pressure.
        '''

        # Enable animations and animation saving in the general sense.
        self.enable_anim_save()

        # Localize nested dictionaries for convenience.
        results = self._config['results options']
        variable = self._config['variable settings']

        # Enable optional features improving test coverage but *NOT* required by
        # the plots and animations enabled below.
        results['enumerate cells'] = True
        results['overlay currents'] = True

        # Enable all plots.
        results['Vmem 2D']['plot Vmem'] = True
        results['Ca 2D']['plot Ca'] = True
        results['pH 2D']['plot pH'] = True
        results['Efield 2D']['plot Efield'] = True
        results['Currents 2D']['plot Currents'] = True
        results['Pressure 2D']['plot Pressure'] = True
        results['Velocity 2D']['plot Velocity'] = True

        # Enable all animations.
        results['Vmem Ani']['animate Vmem'] = True
        results['Ca Ani']['animate Ca2+'] = True
        results['pH Ani']['animate pH'] = True
        results['Vmem GJ Ani']['animate Vmem with gj'] = True
        results['Current Ani']['animate current'] = True
        results['Membrane Ani']['animate membrane'] = True
        results['Efield Ani']['animate Efield'] = True
        results['Velocity Ani']['animate Velocity'] = True
        results['Deformation Ani']['animate Deformation'] = True

        # Enable all features required by these plots and animations.
        self.is_ecm = True
        self.ion_profile = 'animal'
        variable['channel electroosmosis']['turn on'] = True
        variable['deformation']['turn on'] = True
        variable['fluid flow']['include fluid flow'] = True
        # variable['pressures']['include osmotic pressure'] = True # FIXME check this!
