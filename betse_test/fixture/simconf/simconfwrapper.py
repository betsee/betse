#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Test-specific simulation configuration classes wrapping low-level dictionaries
both serialized to and deserialized from on-disk YAML-formatted files.
'''

#FIXME: This submodule requires heavy refactoring away from the current
#low-level approach in favour of the new high-level "confabc"-based approach --
#namely, the "YamlABC" base class coupled with the conf_alias() data
#descriptor. Together, this new functionality provides a substantially better
#YAML-to-Python-object-mapping (YPOM) than the ad-hoc and overly verbose
#boilerplate implemented below. Specifically:
#
#* Eliminate all properties defined below.
#* Replace all usage of the low-level "self.p.conf" dictionary with high-level
#  data descriptors defined by "self.p".
#FIXME: Ideally, after implementing the above, use of the
#"SimConfigTestWrapper" wrapper will be able to be replaced everywhere in tests
#by direct use of the "SimConfTestInternal.p" property providing direct access
#to this "Parameters" object, which increasingly provides all test
#functionality. We're not quite there yet -- but we will be, eventually.

# ....................{ IMPORTS                           }....................
# This subclass imports from submodules defined by the main codebase and is
# thus *NOT* safely importable from fixture submodules directly imported by
# "conftest" plugin modules. To defer the importation of this submodule until
# *AFTER* test collection, this submodule is intentionally segregated.
from betse.science.enum.enumconf import IonProfileType, SolverType
from betse.science.enum.enumphase import SimPhaseKind
from betse.science.parameters import Parameters
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.require import phasereqs
from betse.science.phase.require.abc.phasereqset import (
    SimPhaseRequirementsOrNoneTypes)
from betse.science.pipe.export.pipeexps import SimPipesExport
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ SUPERCLASSES                      }....................
class SimConfigTestWrapper(object):
    '''
    **Test-specific simulation configuration wrapper** (i.e., object wrapping a
    Python dictionary or sequence deserialized from a YAML-formatted simulation
    configuration file).

    This wrapper is intended to be instantiated *only* by non-interactive
    automation (e.g., tests, scripts).

    This wrapper superficially wraps this Python container with convenience
    methods safely modifying this container. This wrapper is lower level than
    the high-level :class:`Parameters` simulation configuration, which
    transforms this container into numerous high-level objects rather than
    merely wrapping this container. While the :class:`Parameters` configuration
    is principally used by backend simulation modelling in the
    :mod:`betse.science` package, this wrapper is principally used by frontend
    logic modifying simulation configurations on behalf of either interactive
    users (e.g., the BETSEE GUI) or automated tests.

    Attributes
    ----------
    p : Parameters
        High-level simulation configuration encapsulated by this test wrapper.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self, p: Parameters) -> None:
        '''
        Wrap the passed simulation configuration.

        Parameters
        ----------
        p : Parameters
            Simulation configuration to be wrapped.
        '''

        # Classify all passed parameters.
        self.p = p

    # ..................{ PROPERTIES ~ float                }..................
    @property
    def environment_size(self) -> float:
        '''
        Square dimension in meters of the external environment containing the
        cell cluster for this configuration.

        For simplicity, BETSE constrains the environment to be square in shape.
        This dimension thus defines both the environmental width *and* height.
        '''

        # Coerce the current number to a float for safety.
        return float(self.p.conf['world options']['world size'])


    @environment_size.setter
    @type_check
    def environment_size(self, environment_size: float) -> None:
        '''
        Set the square dimension in meters of the external environment
        containing the cell cluster for this configuration.
        '''

        # Coerce the passed number to a float for safety.
        self.p.conf['world options']['world size'] = float(environment_size)

    # ..................{ MINIMIZERS                        }..................
    #FIXME: Additionally, the three time durations defined by the
    #"self.p.conf['gene regulatory network settings']['sim-grn settings']"
    #dictionary must also be minified. For simplicity, unconditionally minify
    #these settings within this method regardless of whether the "sim-grn"
    #subcommand is being exercised by the current test or fixture.
    def minify(self) -> None:
        '''
        Minimize the space and time costs associated with running the
        simulation configured by this configuration while preserving all
        fundamental configuration features.

        Specifically, this method numerically reduces all configuration options
        pertaining to either world size *or* simulation time to their **minimum
        permissible values** (i.e., the smallest values still preserving
        simulation stability). This method is intended to be called only by
        test automation.
        '''

        # Minify initialization time to exactly three sampled time steps. For
        # safety, permit currently defined options smaller than the minimums
        # defined below to override these minimums.

        # Duration of each time step in seconds.
        self.p.init_time_step = min(self.p.init_time_step, 1.0e-3)

        # Interval to sample these steps at in seconds.
        self.p.init_time_sampling = min(
            self.p.init_time_sampling, self.p.init_time_step)

        # Total simulation time in seconds. The first digit effectively defines
        # the number of sampled time steps, by the above choice of time step.
        self.p.init_time_total = self.p.init_time_step * 3

        # Minify simulation time to the same durations. To ensure that
        # property setter validation compares durations in the expected manner,
        # properties are assigned in order of increasing duration.
        self.p.sim_time_step = min(
            self.p.sim_time_step, self.p.init_time_step)
        self.p.sim_time_sampling = min(
            self.p.sim_time_sampling, self.p.sim_time_step)
        self.p.sim_time_total = self.p.sim_time_step * 3

        #FIXME: Generalize to minify the time steps of *ALL* enabled events.
        # Minify simulation event times to the same durations in a manner
        # preserving the relative time that these events occur with respect to
        # the pre-minified simulation duration in seconds. For simplicity, each
        # such minification preserves the following ratio:
        #
        #     event_time_old       event_time_new
        #     ------------------ = ------------------
        #     sim_time_total_old   sim_time_total_new
        #
        # Since "sim_time_total_new" is given by "self.p.sim_time_total", the
        # desired unknown "event_time_new" is given by:
        #
        #                      self.p.sim_time_total * event_time_old
        #                      ---------------------------------------
        #     event_time_new = sim_time_total_old
        #FIXME: Well, we've eliminated the "event_cut['cut time']" key entirely.
        #Nonetheless, the above commentary still holds for the general case.

        # Pre-minified simulation duration in seconds.
        # sim_time_total_old = self.p.sim_time_total

        # event_cut = self.p.conf['cutting event']
        # event_cut['cut time'] = (
        #     self.p.sim_time_total * event_cut['cut time'] /
        #     sim_time_total_old)

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
        self.p.grid_size = min(self.p.grid_size, 20)

        # Log this minification.
        logs.log_debug(
            'Minifying simulation to init '
            'for %f s at [ %f | %f ] s and sim '
            'for %f s at [ %f | %f ] s...',
            self.p.init_time_step,
            self.p.init_time_sampling,
            self.p.init_time_total,
            self.p.sim_time_step,
            self.p.sim_time_sampling,
            self.p.sim_time_total,
        )

    # ..................{ DISABLERS                         }..................
    #FIXME: The implementation of the following methods is fundamentally unsafe.
    #If the structure of the underlying YAML file changes, these methods could
    #silently fail (e.g., if the "plot while solving" option were renamed to
    #"is plotting during"). To combat this, all attempts to directly modify the
    #"self.p.conf" dictionary below *MUST* instead defer to a newly defined
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

        self.p.anim.is_after_sim = False
        self.p.anim.is_while_sim = False
        self.p.plot.is_after_sim = False


    def disable_interaction(self) -> None:
        '''
        Disable all simulation configuration options either requiring
        interactive user input _or_ displaying graphical output intended for
        interactive user consumption (e.g., plots, animations).

        This method is intended to be called by non-interactive automation
        (e.g., tests, scripts) expecting simulations to behave silently.
        '''

        self.p.anim.is_after_sim_show = False
        self.p.anim.is_while_sim_show = False
        self.p.plot.is_after_sim_show = False

    # ..................{ ENABLERS                          }..................
    def enable_networks(self) -> None:
        '''
        Enable both biochemical reaction and gene regulatory networks.
        '''

        self.p.conf['gene regulatory network settings'][
            'gene regulatory network simulated'] = True

    # ..................{ ENABLERS ~ export                 }..................
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
            (e.g., ``ffmpeg``, ``imagemagick``).
        filetype : str
            Filetype of videos to encode with this writer (e.g., ``mkv``).
        '''

        # Enable animations and animation saving in the general sense.
        self.enable_visuals_save()

        # Localize nested dictionaries for convenience.
        video = self.p.conf['results options']['save']['animations']['video']

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

        self.p.anim.is_while_sim_save = True
        self.p.anim.is_after_sim_save = True
        self.p.plot.is_after_sim_save = True

    # ..................{ ENABLERS ~ solver : fast          }..................
    def _enable_solver_fast(self) -> None:
        '''
        Enable the equivalent circuit-based BETSE solver *and* disable all
        simulation features unsupported by this solver.
        '''

        self.p.solver_type = SolverType.FAST


    def enable_solver_fast_exports(self) -> None:
        '''
        Enable all possible exports (e.g., CSVs, plots, animations) supported
        by the equivalent circuit-based BETSE solver excluding those requiring
        extracellular spaces, all features required by these exports, and any
        additional features trivially enabled *without* increasing time or
        space complexity.
        '''

        # Enable all simulation features, including the full BETSE solver but
        # excluding extracellular spaces.
        self._enable_solver_fast_features()

        # Disable extracellular spaces.
        self.p.is_ecm = False

        # Enable all possible exports excluding those requiring the full
        # solver.
        self._enable_exports(requirements_omit=phasereqs.SOLVER_FULL)


    def _enable_solver_fast_features(self) -> None:
        '''
        Enable all simulation features required by all exports (e.g., CSVs,
        plots, animations) supported by the equivalent circuit-based BETSE
        solver, excluding extracellular spaces.

        This method additionally enables optional settings improving test
        coverage but *not* explicitly required by these exports. Specifically,
        this method enables:

        * The full BETSE solver.
        * Saving of all visual exports.
        '''

        # Enable the equivalent circuit-based solver.
        self._enable_solver_fast()

        # Enable saving of these exports.
        self.enable_visuals_save()

    # ..................{ ENABLERS ~ solver : full          }..................
    def _enable_solver_full(self) -> None:
        '''
        Enable the complete BETSE solver.
        '''

        self.p.solver_type = SolverType.FULL


    def enable_solver_full_vg_ions(self) -> None:
        '''
        Enable all voltage-gated ion channels (e.g., sodium, potassium) *and*
        all features required by these channels.

        This method is intended to be called by non-interactive test suites
        exercising these channels. Specifically, this method enables:

        * The full BETSE solver.
        * The extracellular matrix (ECM).
        * The mammalian ion profile (i.e., ``animal``), enabling all ions.
        * The intervention increasing the permeability of all cell membranes to
          sodium (Na+).
        * The voltage-gated sodium (Na+) channel ``Nav1p2``, corresponding to
          the adult human brain.
        * The voltage-gated potassium (K+) channel ``K_Slow``.
        * Decreased time step and sampling rates, ensuring simulation
          stability.
        * Increased duration and cell count, exposing simulation instabilities.

        For efficiency, this method disables all visuals -- including both in-
        and post-simulation animations and plots.
        '''

        # For efficiency, disable all visuals.
        self.disable_visuals()

        # Enable all features required by these channels.
        self._enable_solver_full()
        self.p.is_ecm = True
        self.p.ion_profile = IonProfileType.MAMMAL

        # For stability, decrease both the time step and sampling rates.
        self.p.sim_time_step     = 1e-4
        self.p.sim_time_sampling = 1e-3

        # For completeness, increase both the duration and cell count.
        self.p.sim_time_total = 50e-3
        self.environment_size = 250e-6

        # Enable the intervention increasing sodium membrane permeability.
        # Although the current default values for this intervention track those
        # defined below fairly closely, the latter are nonetheless explicitly
        # defined below to avoid issues when the former inevitably change.
        sodium_membrane_permeability = self.p.conf['change Na mem']
        sodium_membrane_permeability['event happens'] = True
        sodium_membrane_permeability['change rate']   =  1.0e-3
        sodium_membrane_permeability['change start']  =  5.0e-3
        sodium_membrane_permeability['change finish'] = 30.0e-3
        sodium_membrane_permeability['apply to'] = ['Spot',]

        #FIXME: Refactor this to use the new networks formalism.
        # # Enable the voltage-gated sodium (Na+) channel Nav1p2.
        # voltage_gated_sodium_channel = self.p.conf['voltage gated Na+']
        # voltage_gated_sodium_channel['turn on'] = True
        # voltage_gated_sodium_channel['channel type'] = ['Nav1p2',]
        # # voltage_gated_sodium_channel['max value'] = 5.0e-6
        # voltage_gated_sodium_channel['apply to'] = ['base',]
        #
        # # Enable the voltage-gated potassium (K+) channel K_Slow.
        # voltage_gated_potassium_channel = self.p.conf['voltage gated K+']
        # voltage_gated_potassium_channel['turn on'] = True
        # voltage_gated_potassium_channel['channel type'] = ['K_Slow',]
        # # voltage_gated_potassium_channel['max value'] = 5.0e-7
        # voltage_gated_potassium_channel['apply to'] = ['base',]

    # ..................{ ENABLERS ~ solver : full : exports }..................
    def enable_solver_full_exports_ecm(self) -> None:
        '''
        Enable all possible exports (e.g., CSVs, plots, animations) supported
        by the full BETSE solver including those requiring extracellular spaces,
        all features required by these exports, and any additional features
        trivially enabled *without* increasing time or space complexity.
        '''

        # Enable all simulation features, including the full BETSE solver but
        # excluding extracellular spaces.
        self._enable_solver_full_features()

        # Enable extracellular spaces.
        self.p.is_ecm = True

        # Enable all possible exports.
        self._enable_exports()


    def enable_solver_full_exports_noecm(self) -> None:
        '''
        Enable all possible exports (e.g., CSVs, plots, animations) supported
        by the full BETSE solver excluding those requiring extracellular spaces,
        all features required by these exports, and any additional features
        trivially enabled *without* increasing time or space complexity.
        '''

        # Enable all simulation features, including the full BETSE solver but
        # excluding extracellular spaces.
        self._enable_solver_full_features()

        # Disable extracellular spaces.
        self.p.is_ecm = False

        # Enable all possible exports excluding those requiring extracellular
        # spaces.
        self._enable_exports(requirements_omit=phasereqs.ECM)

    # ..................{ PRIVATE ~ enablers                }..................
    @type_check
    def _enable_exports(
        self,

        # Optional parameters.
        requirements_omit: SimPhaseRequirementsOrNoneTypes = None,
    ) -> None:
        '''
        Enable all possible exports (e.g., CSVs, plots, animations) excluding
        those requiring one or more of the passed simulation phase requirements
        (e.g., the full BETSE solver, extracellular spaces).

        Exports already enabled by the current simulation configuration are
        gracefully preserved as is rather than erroneously readded to their
        respective pipelines.

        Parameters
        ----------
        requirements_omit : SimPhaseRequirementsOrNoneTypes
            Immutable set of all simulation phase requirements such that
            exports requiring one or more requirements in this set are *not*
            enabled by this method. Defaults to ``None``, in which case all
            possible exports are unconditionally enabled.
        '''

        # Defer heavyweight imports.
        from betse.science.config.export.confexpabc import SimConfExportABC
        from betse.util.type.obj import objtest

        # Log this action.
        logs.log_debug('Analyzing pipeline exporters...')

        # Simulation phase encapsulating this configuration. Since this phase
        # *ONLY* serves as a thin wrapper around this configuration, the type
        # of phase is irrelevant; ergo, this type arbitrarily defaults to the
        # first phase type: seed.
        phase = SimPhase(kind=SimPhaseKind.SEED, p=self.p)

        # Default the set of requirements to omit to the empty set.
        if requirements_omit is None:
            requirements_omit = phasereqs.NONE

        # For each export pipeline...
        for pipe_export in SimPipesExport().PIPES_EXPORT:
            # Log this pipeline.
            logs.log_debug(
                'Analyzing pipeline "%s" exporters...', pipe_export.name)

            # Sequence of all export subconfigurations for this pipeline.
            pipe_exporters_conf = pipe_export.iter_runners_conf(phase)

            # Remove all export subconfigurations from this sequence, enabling
            # new test-specific subconfigurations to be added to this sequence
            # below without needless concern over conflicts and redundancy.
            pipe_exporters_conf.clear()

            # For each possible export supported by this pipeline (including
            # those currently disabled)...
            for pipe_exporter_metadata in pipe_export.iter_runners_metadata():
                # Log this iteration.
                logs.log_debug(
                    'Analyzing pipeline "%s" exporter "%s"...',
                    pipe_export.name, pipe_exporter_metadata.kind)

                # If this export requires one or more simulation features
                # omitted by the caller...
                if pipe_exporter_metadata.requirements.isintersection(
                    requirements_omit):
                    # Log this exclusion.
                    logs.log_debug(
                        'Excluding pipeline "%s" exporter "%s", '
                        'due to unsatisfied test requirements...',
                        pipe_export.name, pipe_exporter_metadata.kind)

                    # Continue to the next export.
                    continue
                # Else, this export requires no omitted simulation features.

                # Log this inclusion.
                logs.log_debug(
                    'Including pipeline "%s" exporter "%s"...',
                    pipe_export.name, pipe_exporter_metadata.kind)

                # Create and append a new default subconfiguration of this
                # export to the sequence of these subconfigurations.
                pipe_exporter_conf = pipe_exporters_conf.append_default()

                # If this is *NOT* actually an export subconfiguration, raise
                # an exception.
                objtest.die_unless_instance(
                    obj=pipe_exporter_conf, cls=SimConfExportABC)

                # Copy across the type of this export subconfiguration.
                pipe_exporter_conf.kind = pipe_exporter_metadata.kind

    # ..................{ PRIVATE ~ enablers : solver       }..................
    def _enable_solver_full_features(self) -> None:
        '''
        Enable all simulation features required by all exports (e.g., CSVs,
        plots, animations) -- including the full BETSE solver but excluding
        extracellular spaces.

        This method additionally enables optional settings improving test
        coverage but *not* explicitly required by these exports. Specifically,
        this method enables:

        * The full BETSE solver.
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

        # Enable the full solver.
        self._enable_solver_full()

        # Enable saving of these exports.
        self.enable_visuals_save()

        # Localize nested dictionaries for convenience.
        results = self.p.conf['results options']
        variable = self.p.conf['variable settings']

        # Enable all simulation features required by these exports.
        self.p.ion_profile = IonProfileType.MAMMAL
        self.p.conf['apply pressure']['event happens'] = True
        # variable['channel electroosmosis']['turn on'] = True  # This feature has been removed
        variable['deformation']['turn on'] = True
        # variable['fluid flow']['include fluid flow'] = True
        variable['pressures']['include osmotic pressure'] = True
        self.p.cell_polarizability = 1e-4

        # Enable all optional settings supported by these exports.
        results['visuals']['cell indices']['show'] = True
        results['overlay currents'] = True
