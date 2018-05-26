#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from betse.exceptions import (
    BetseFileException, BetseSimException, BetseSimConfException)
from betse.lib.pickle import pickles
from betse.science import filehandling as fh
from betse.science.cells import Cells
from betse.science.chemistry.gene import MasterOfGenes
from betse.science.config.confenum import GrnUnpicklePhaseType
from betse.science.export import exppipe
from betse.science.parameters import Parameters
from betse.science.sim import Simulator
from betse.science.phase import phasecallbacks
from betse.science.phase.phasecallbacks import SimCallbacksABCOrNoneTypes
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.phaseenum import SimPhaseKind
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.type.decorator.decorators import deprecated
from betse.util.type.decorator.decprof import log_time_seconds
from betse.util.type.types import type_check

# ....................{ CLASSES                           }....................
class SimRunner(object):
    '''
    **Simulation runner** (i.e., high-level object encapsulating the running of
    simulation phases as corresponding public methods commonly referred to as
    simulation subcommands).

    This runner provides high-level methods for initializing, running, and
    plotting simulations specified by the YAML-formatted simulation
    configuration file with which this runner is instantiated. Each instance of
    this runner simulates exactly one simulation.

    Attributes
    ----------
    _callbacks : SimCallbacksABC
        Caller-defined object whose methods are periodically called during each
        simulation subcommand (e.g., to notify this caller of phase progress).
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        p: Parameters,

        # Optional parameters.
        callbacks: SimCallbacksABCOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this simulation runner.

        Attributes
        ----------
        p : Parameters
            Simulation configuration to be run. If this configuration is in the
            unloaded state (i.e., has *not* yet been loaded into memory from a
            YAML-formatted file on disk), an exception is raised.
        callbacks : SimCallbacksABCOrNoneTypes
            Caller-defined object whose methods are periodically called during
            each simulation subcommand (e.g., :meth:`SimRunner.seed`). Defaults
            to ``None``, in which case this defaults to a placeholder object
            whose methods all silently reduce to noops.

        Raises
        ----------
        BetseYamlException
            If this simulation configuration is in the unloaded state.
        '''

        # If this configuration is unloaded, raise an exception.
        p.die_unless_loaded()

        # Default all unpassed parameters to sane defaults.
        if callbacks is None:
            callbacks = phasecallbacks.make_default()

        # Classify all passed parameters *AFTER* defaulting these parameters.
        self._callbacks = callbacks
        self._p = p

    # ..................{ RUNNERS                            }..................
    @log_time_seconds(noun='seed')
    def seed(self) -> SimPhase:
        '''
        Seed this simulation with a new cell cluster and cache this cluster to
        an output file, specified by the current configuration file.

        This method *must* be called prior to the :meth:`init` and
        :meth:`plot_seed` methods, which consume this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Cuumulative number of times that each call of the this subcommand
        # calls the SimCallbacksABC.progressed() callback.
        #
        # This magic number *MUST* be manually synchronized with the
        # implementation of both this method and methods transitively called by
        # this method. Failure to do so *WILL* raise exceptions. Sadly, there
        # exists no reasonable means of automating this synchronization.
        SEED_PROGRESS_TOTAL = 6

        # Log this attempt.
        logs.log_info('Seeding simulation...')

        # Notify the caller of the range of work performed by this subcommand.
        self._callbacks.progress_ranged(progress_max=SEED_PROGRESS_TOTAL)

        # Simulation phase.
        phase = SimPhase(
            kind=SimPhaseKind.SEED, p=self._p, callbacks=self._callbacks)

        # Create the pseudo-randomized cell cluster.
        self._callbacks.progressed_next()
        phase.cells.make_world(phase)

        # Initialize core simulation data structures.
        self._callbacks.progressed_next()
        phase.sim.init_core(phase)

        # Define the tissue and boundary profiles for plotting.
        self._callbacks.progressed_next()
        phase.dyna.tissueProfiles(phase.sim, phase.cells, self._p)

        # Redo gap junctions to isolate different tissue types.
        self._callbacks.progressed_next()
        phase.cells.redo_gj(phase.dyna, self._p)

        # Create a Laplacian and solver for discrete transfers on closed,
        # irregular cell network.
        self._callbacks.progressed_next()
        phase.cells.graphLaplacian(self._p)

        #FIXME: Would shifting this logic into the cells.graphLaplacian() method
        #called above be feasible? If not, no worries! (Granular lunar sunsets!)
        if not self._p.td_deform:  # if time-dependent deformation is not required
            phase.cells.lapGJ = None
            phase.cells.lapGJ_P = None  # null out the non-inverse matrices -- we don't need them

        # Create accessory matrices depending on user requirements.
        if self._p.deformation:
            phase.cells.deform_tools(self._p)

        # if self._p.sim_eosmosis is True:
        #     phase.cells.eosmo_tools(self._p)

        # Finish up.
        self._callbacks.progressed_next()
        phase.cells.save_cluster(self._p)

        # Log the completion of this phase.
        logs.log_info('Cell cluster creation complete!')
        phase.sim.sim_info_report(phase.cells, self._p)

        # Return this phase.
        return phase


    @log_time_seconds(noun='initialization')
    def init(self) -> SimPhase:
        '''
        Initialize this simulation with the cell cluster seeded by a prior call
        to the :meth:`seed` method and cache this initialization to an output
        file, specified by the current configuration file.

        This method *must* be called prior to the :meth:`sim` and
        :meth:`plot_init` methods, which consume this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Log this attempt.
        logs.log_info('Initializing simulation...')

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # High-level simulation objects.
        cells = Cells(self._p)  # create an instance of world
        sim = Simulator(self._p)

        #FIXME: The Parameters.__init__() method should *REQUIRE* that a time
        #profile type be passed. The current approach leaves critical attributes
        #undefined in the event that the optional Parameters.set_time_profile()
        #method is left uncalled, which is pretty unacceptable.
        #FIXME: Actually, no. All logic performed by the set_time_profile()
        #method should be shifted into the SimPhase.__init__() method. See a
        #FIXME comment preceding the set_time_profile() method for details.
        self._p.set_time_profile(phase_kind)  # force the time profile to be initialize
        self._p.run_sim = False # let the simulator know we're just running an initialization
        sim.run_sim = False

        #FIXME: This if conditional is repeated verbatim twice below. Generalize
        #into a new _load_cells() method containing this if conditional and
        #returning the loaded "Cells" instance; then, call this method both here
        #and everywhere this repeated logic appears below. Starbust dragons!
        if files.is_file(cells.savedWorld):
            cells, p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')

            # check to ensure compatibility between original and present sim files:
            self._die_if_seed_differs(p_old, self._p)

        else:
            logs.log_warning("Ooops! No such cell cluster file found to load!")

            if self._p.autoInit:
                logs.log_info(
                    'Automatically seeding cell cluster from config file settings...')
                self.seed()  # create an instance of world
                logs.log_info(
                    'Now using cell cluster to run initialization.')
                cells,_ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

            else:
                raise BetseSimException(
                    "Run terminated due to missing seed.\n"
                    "Please run 'betse seed' to try again.")

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        sim.init_core(phase)

        # Run this simulation phase.
        sim.sim_info_report(cells, self._p)
        sim.run_sim_core(phase)

        # Return this phase.
        return phase


    @log_time_seconds(noun='simulation')
    def sim(self) -> SimPhase:
        '''
        Simulate this simulation with the cell cluster initialized by a prior
        call to the :meth:`init` method and cache this simulation to an output
        file, specified by the current configuration file.

        This method *must* be called prior to the :meth:`:meth:`plot_sim`
        method, which consumes this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Log this attempt.
        logs.log_info('Running simulation...')

        # Simulation phase type.
        phase_kind = SimPhaseKind.SIM

        # High-level simulation objects.
        sim = Simulator(p=self._p)

        #FIXME: See above for pertinent commentary. Mists of time, unpart!
        self._p.set_time_profile(phase_kind)  # force the time profile to be initialize
        self._p.run_sim = True    # set on the fly a boolean to let simulator know we're running a full simulation

        if files.is_file(sim.savedInit):
            # Load the initialization from cache.
            sim, cells, p_old = fh.loadSim(sim.savedInit)

            # Ensure compatibility between original and present config files.
            self._die_if_seed_differs(p_old, self._p)
        else:
            logs.log_warning(
                "No initialization file found to run this simulation!")

            if self._p.autoInit:
                logs.log_info("Automatically running initialization...")
                self.init()
                logs.log_info('Now using initialization to run simulation.')

                # Load the initialization from cache.
                sim, cells, _ = fh.loadSim(sim.savedInit)

            else:
                raise BetseSimException(
                    'Simulation terminated due to missing initialization. '
                    'Please run an initialization and try again.')

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        # Run and save the simulation to the cache.
        sim.sim_info_report(cells, self._p)
        sim.run_sim_core(phase)

        # Return this phase.
        return phase


    @log_time_seconds(noun='network')
    def sim_grn(self) -> SimPhase:
        '''
        Initialize and simulate a pure gene regulatory network (GRN) *without*
        bioelectrics with the cell cluster seeded by a prior call to the
        :meth:`seed` method and cache this initialization and simulation to
        output files, specified by the current configuration file.

        This method *must* be called prior to the :meth:`plot_grn` method, which
        consumes this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Simulation phase objects, defaulting to undefined initially.
        phase = None

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        #FIXME: Non-ideal. Ideally, these objects would both default to "None".
        #Sadly, they're currently required for trivial access to filename
        #variables (e.g., "cells.savedWorld", "sim.savedInit"). After shifting
        #all such variables into the "Parameters" class, uncomment the
        #following two lines and remove the corresponding two lines below.
        # cells = None
        # sim   = None

        # Simulation configuration.
        cells = Cells(self._p)
        sim = Simulator(self._p)

        # Log this simulation.
        logs.log_info(
            'Running gene regulatory network "%s" '
            'defined in config file "%s"...',
            pathnames.get_basename(self._p.grn_config_filename),
            self._p.conf_basename)

        #FIXME: See above for pertinent commentary. Tendrils of wisdom, uncoil!
        self._p.set_time_profile(phase_kind)  # force the time profile to be initialize
        self._p.run_sim = False

        # If networking an uninitialized, unsimulated cell cluster...
        if self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.SEED:
            if files.is_file(cells.savedWorld):
                cells, _ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
                logs.log_info('Running gene regulatory network on betse seed...')

                # Simulation phase.
                phase = SimPhase(
                    kind=phase_kind,
                    cells=cells,
                    p=self._p,
                    sim=sim,
                    callbacks=self._callbacks,
                )

                # Initialize core simulation data structures.
                sim.init_core(phase)
                sim.init_dynamics(phase)

                # Initialize other aspects required for piggyback of GRN on the
                # sim object.
                sim.time = []
                sim.vm = -50e-3 * np.ones(sim.mdl)

                # Initialize key fields of simulator required to interface
                # (dummy init).
                sim.rho_pump = 1.0
                sim.rho_channel = 1.0
                sim.conc_J_x = np.zeros(sim.edl)
                sim.conc_J_y = np.zeros(sim.edl)
                sim.J_env_x = np.zeros(sim.edl)
                sim.J_env_y = np.zeros(sim.edl)
                sim.u_env_x = np.zeros(sim.edl)
                sim.u_env_y = np.zeros(sim.edl)

            else:
                logs.log_warning(
                    'Ooops! No such cell cluster file found to load!')

                if self._p.autoInit:
                    logs.log_info(
                        'Automatically seeding cell cluster from config file settings...')
                    self.seed()  # create an instance of world
                    logs.log_info(
                        'Now using cell cluster to run initialization.')
                    cells, _ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache
                else:
                    raise BetseSimException(
                        'Run terminated due to missing seed. '
                        'Please run "betse seed" to try again.')
        # Else if networking an initialized but unsimulated cell cluster...
        elif self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.INIT:
            if files.is_file(sim.savedInit):
                logs.log_info('Running gene regulatory network on betse init...')
                sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
            else:
                logs.log_warning(
                    "No initialization file found to run the GRN simulation!")

                if self._p.autoInit:
                    logs.log_info("Automatically running initialization...")
                    self.init()
                    logs.log_info('Now using initialization to run simulation.')
                    sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
                else:
                    raise BetseSimException(
                        'Simulation terminated due to missing core initialization. '
                        'Please run a betse initialization and try again.')
        # Else if networking an initialized, simulated cell cluster...
        elif self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.SIM:
            if files.is_file(sim.savedSim):
                logs.log_info('Running gene regulatory network on betse sim...')
                sim, cells, _ = fh.loadSim(sim.savedSim)  # load the initialization from cache
            else:
                logs.log_warning(
                    "No simulation file found to run the GRN simulation!")
                raise BetseSimException(
                    'Simulation terminated due to missing core simulation. '
                    'Please run a betse simulation and try again.')
        # Else, this type of networking is unrecognized. Raise an exception.
        else:
            raise BetseSimConfException(
                'Gene regulatory network (GRN) unpickle simulation phase '
                '"{}" unrecognized.'.format(self._p.grn_unpickle_phase_type))

        # If *NOT* defined above, define this simulation phase.
        if phase is None:
            phase = SimPhase(
                kind=phase_kind,
                cells=cells,
                p=self._p,
                sim=sim,
                callbacks=self._callbacks,
            )

            # Reinitialize all profiles.
            phase.dyna.tissueProfiles(sim, cells, self._p)

        # If *NOT* restarting from a prior GRN run, start a new GRN.
        if self._p.grn_unpickle_filename is None:
            # Log this start.
            logs.log_info("Initializing the gene regulatory network...")

            # Create and initialize an instance of master of metabolism.
            MoG = MasterOfGenes(self._p)
            MoG.read_gene_config(sim, cells, self._p)
        # Else, restart from a prior GRN run.
        else:
            # Log this restart.
            logs.log_info(
                'Reinitializing the gene regulatory network from "%s"...',
                pathnames.get_basename(self._p.grn_unpickle_filename))

            # If this file does *NOT* exist, raise an exception.
            if not files.is_file(self._p.grn_unpickle_filename):
                raise BetseSimException(
                    'Gene regulatory network unloadable '
                    'from file not found: {}'.format(
                        self._p.grn_unpickle_filename))

            # Unpickle this file into a high-level "MasterOfGenes" object.
            MoG, _, _ = pickles.load(self._p.grn_unpickle_filename)

            # If running on a sim with a cut event, perform this cut...
            if (
                phase.dyna.event_cut is not None and
                phase.dyna.event_cut.is_fired
            ):
                # Log this cutting.
                logs.log_info(
                    'A cutting event has been run, '
                    'so the GRN object needs to be modified...')

                sim_old = Simulator(self._p)

                # If no prior initialization exists, raise an exception.
                if not files.is_file(sim_old.savedInit):
                    logs.log_warning(
                        'This situation is complex '
                        'due to a cutting event being run. '
                        'Please have a corresponding init file '
                        'to run the GRN simulation!')

                    raise BetseSimException(
                        'Simulation terminated due to missing core init. '
                        'Please alter GRN settings and try again.')
                # Else, a prior initialization exists.

                # Log this initialization.
                logs.log_info(
                    'Loading betse init from cache '
                    'for reference to original cells...')

                # Load the initialization from cache.
                init, cells_old, _ = fh.loadSim(sim_old.savedInit)

                #FIXME: This phase object would ideally be pickled to and
                #from the "sim_old.savedInit" file loaded above, in which
                #case this local variable would be safely removable. Flagon!

                # Original simulation phase. To avoid caller confusion, the
                # optional "callbacks" parameter is intentionally *NOT* passed.
                phase_old = SimPhase(
                    kind=phase_kind,
                    cells=cells_old,
                    p=self._p,
                    sim=sim_old,
                )

                # Initialize all tissue profiles on original cells.
                phase_old.dyna.tissueProfiles(init, cells_old, self._p)

                for cut_profile_name in (
                    phase_old.dyna.event_cut.profile_names):
                    logs.log_info(
                        'Cutting cell cluster via cut profile "%s"...',
                        cut_profile_name)

                    # Object picking the cells removed by this cut profile.
                    tissue_picker = phase_old.dyna.cut_name_to_profile[
                        cut_profile_name].picker

                    # One-dimensional Numpy arrays of the indices of all
                    # cells and cell membranes to be removed.
                    target_inds_cell, target_inds_mem = (
                        tissue_picker.pick_cells_and_mems(
                            cells=cells_old, p=self._p))

                    MoG.core.mod_after_cut_event(
                        target_inds_cell,
                        target_inds_mem,
                        sim,
                        cells,
                        self._p,
                    )

                    logs.log_info(
                        'Redefining dynamic dictionaries '
                        'to point to the new sim...')
                    MoG.core.redefine_dynamic_dics(sim, cells, self._p)

                    logs.log_info(
                        'Reinitializing the gene regulatory network '
                        'for simulation...')
                    MoG.reinitialize(sim, cells, self._p)

            #FIXME: Would you mind documenting this a bit, beautiful! Lovalo!
            if self._p.use_microtubules:
                sim.mtubes.mtubes_x = MoG.mtubes_x_time[-1]
                sim.mtubes.mtubes_y = MoG.mtubes_y_time[-1]

                sim.mtubes.uxmt, sim.mtubes.uymt = sim.mtubes.mtubes_to_cell(
                    cells, self._p)

        logs.log_info("Running gene regulatory network test simulation...")
        MoG.run_core_sim(sim, cells, self._p)

        # Return this phase.
        return phase

    # ..................{ PLOTTERS                           }..................
    @log_time_seconds(noun='seed', verb='exported')
    def plot_seed(self) -> SimPhase:
        '''
        Visualize the cell cluster seed by a prior call to the :meth:`seed`
        method and export the resulting plots and animations to various output
        files, specified by the current configuration file.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        logs.log_info(
            'Plotting cell cluster with configuration file "%s".',
            self._p.conf_basename)

        # High-level simulation objects.
        cells = Cells(self._p)
        sim = Simulator(self._p)

        if files.is_file(cells.savedWorld):
            cells, _ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')
        else:
            raise BetseSimException(
                "Ooops! No such cell cluster file found to load!")

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=SimPhaseKind.SEED,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        sim.init_core(phase)
        phase.dyna.tissueProfiles(sim, cells, self._p)

        #FIXME: Refactor into a seed-specific plot pipeline. Dreaming androids!
        if self._p.autosave:
            savedImg = pathnames.join(self._p.init_export_dirname, 'fig_')

        # if self._p.plot_cell_cluster:
        # fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(
        #     self._p, phase.dyna, cells, clrmap=self._p.background_cm)

        if self._p.autosave:
            savename10 = savedImg + 'cluster_mosaic' + '.png'
            plt.savefig(savename10, format='png', transparent=True)

        if self._p.plot.is_after_sim_show:
            plt.show(block=False)

        if self._p.is_ecm:  # and self._p.plot_cluster_mask:
            plt.figure()
            ax99 = plt.subplot(111)
            plt.imshow(
                np.log10(sim.D_env_weight.reshape(cells.X.shape)),
                origin='lower',
                extent=[self._p.um * cells.xmin, self._p.um * cells.xmax, self._p.um * cells.ymin, self._p.um * cells.ymax],
                cmap=self._p.background_cm,
            )
            plt.colorbar()

            cell_edges_flat = self._p.um * cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat, colors='k')
            coll.set_alpha(1.0)
            ax99.add_collection(coll)

            plt.title('Logarithm of Environmental Diffusion Weight Matrix')

            if self._p.autosave:
                savename10 = savedImg + 'env_diffusion_weights' + '.png'
                plt.savefig(savename10, format='png', transparent=True)

            if self._p.plot.is_after_sim_show:
                plt.show(block=False)

            plt.figure()
            plt.imshow(
                cells.maskM,
                origin='lower',
                extent=[self._p.um * cells.xmin, self._p.um * cells.xmax, self._p.um * cells.ymin, self._p.um * cells.ymax],
                cmap=self._p.background_cm,
            )
            plt.colorbar()
            plt.title('Cluster Masking Matrix')

            if self._p.autosave:
                savename = savedImg + 'cluster_mask' + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if self._p.plot.is_after_sim_show:
                plt.show(block=False)

        # Plot gap junctions.
        # if self._p.plot_cell_connectivity:
        plt.figure()
        ax_x = plt.subplot(111)

        if self._p.showCells:
            base_points = np.multiply(cells.cell_verts, self._p.um)
            col_cells = PolyCollection(base_points, facecolors='k', edgecolors='none')
            col_cells.set_alpha(0.3)
            ax_x.add_collection(col_cells)

        con_segs = cells.nn_edges
        connects = self._p.um * np.asarray(con_segs)
        collection = LineCollection(connects, linewidths=1.0, color='b')
        ax_x.add_collection(collection)
        plt.axis('equal')
        plt.axis([cells.xmin * self._p.um, cells.xmax * self._p.um, cells.ymin * self._p.um, cells.ymax * self._p.um])

        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Cell Connectivity Network')

        if self._p.autosave is True:
            savename10 = savedImg + 'gj_connectivity_network' + '.png'
            plt.savefig(savename10, format='png', transparent=True)

        if self._p.turn_all_plots_off is False:
            plt.show(block=False)

        else:
            logs.log_info(
                'Plots exported to init results folder '
                'defined in configuration file "%s".',
                self._p.conf_basename)

        # Return this phase.
        return phase


    @log_time_seconds(noun='initialization', verb='exported')
    def plot_init(self) -> SimPhase:
        '''
        Visualize the cell cluster initialized by a prior call to the
        :meth:`init` method and export the resulting plots and animations to
        various output files, specified by the current configuration file.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Log this plotting attempt.
        logs.log_info(
            'Plotting initialization with configuration "%s"...',
            self._p.conf_basename)

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # High-level simulation objects.
        sim = Simulator(p=self._p)

        #FIXME: See above for pertinent commentary. Mists of time, unpart!
        self._p.set_time_profile(phase_kind)  # force the time profile to be initialize

        #FIXME: Bizarre logic. We create a "Simulator" instance above only to
        #test whether a single file exists and, if so, replace that instance
        #with a pickled "Simulator" instance unpickled from that file. Let's cut
        #out the API middleman, as it were, by obtaining the pathname for this
        #file from the "Parameters" instance instead and then removing the above
        #instantiation of "sim = Simulator(self._p)". Note when doing so that similar
        #behaviour has been duplicated across this submodule. Cheerful cherries!

        if files.is_file(sim.savedInit):
            sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        else:
            raise BetseSimException(
                "Ooops! No such initialization file found to plot!")

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        #FIXME: This... isn't the best. Ideally, the phase.dyna.tissueProfiles()
        #method would *ALWAYS* be implicitly called by the SimPhase.__init__()
        #method. Unfortunately, the non-trivial complexity of cell cluster
        #initialization requires we do so manually for now. Sad sandlion frowns!
        phase.dyna.tissueProfiles(sim, cells, self._p)

        # Display and/or save all initialization exports (e.g., animations).
        exppipe.pipeline(phase)

        #FIXME: All of the following crash if image saving is not turned on, but
        #due to whatever way this is set up, it's not possible to readily fix
        #it. Grrrrrr.....
        #FIXME: Why not just stop each block from happening if image saving is
        #off? Easy peasy!

        #FIXME: Reduce duplication. The following logic is effectively a copy of
        #similar logic in the plot_sim() method.
        #FIXME: Split each of the following blocks performing both plotting and
        #animating into their appropriate plotpipe.pipeline() or
        #animpipe.pipeline() functions, which should resolve the above concerns.

        # run the molecules plots:
        if self._p.molecules_enabled and sim.molecules is not None:
            # reinit settings for plots, in case they've changed:
            sim.molecules.core.plot_init(self._p.network_config, self._p)
            sim.molecules.core.init_saving(cells, self._p, plot_type='init')
            sim.molecules.core.export_all_data(sim, cells, self._p, message='auxiliary molecules')
            sim.molecules.core.plot(sim, cells, self._p, message='auxiliary molecules')
            sim.molecules.core.anim(phase=phase, message='auxiliary molecules')

        if self._p.grn_enabled and sim.grn is not None:
            # reinitialize the plot settings:
            sim.grn.core.plot_init(self._p.grn.conf, self._p)
            sim.grn.core.init_saving(
                cells, self._p, plot_type='init', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, self._p, message='GRN molecules')
            sim.grn.core.plot(sim, cells, self._p, message='GRN molecules')
            sim.grn.core.anim(phase=phase, message='GRN molecules')

        #FIXME: Integrate into the plot pipeline.
        if self._p.Ca_dyn and self._p.ions_dict['Ca'] == 1:
            sim.endo_retic.init_saving(
                cells, self._p, plot_type='init', nested_folder_name='ER')
            sim.endo_retic.plot_er(sim, cells, self._p)

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if self._p.plot.is_after_sim_show:
            plt.show()

        # Return this phase.
        return phase


    @log_time_seconds(noun='simulation', verb='exported')
    def plot_sim(self) -> SimPhase:
        '''
        Visualize the cell cluster simulated by a prior call to the :meth:`sim`
        method and export the resulting plots and animations to various output
        files, specified by the current configuration file.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Log this plotting attempt.
        logs.log_info(
            'Plotting simulation with configuration "%s"...',
            self._p.conf_basename)

        # Simulation phase type.
        phase_kind = SimPhaseKind.SIM

        # High-level simulation objects.
        sim = Simulator(p=self._p)

        #FIXME: See above for pertinent commentary. Mists of time, unpart!
        self._p.set_time_profile(phase_kind)  # force the time profile to be simulation

        # If this simulation has yet to be run, fail.
        if not files.is_file(sim.savedSim):
            raise BetseFileException(
                'Simulation cache file "{}" not found to plot '
                '(e.g., due to no simulation having been run).'.format(
                    sim.savedSim))

        # Load the simulation from the cache.
        sim, cells, _ = fh.loadSim(sim.savedSim)

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        #FIXME: This... isn't the best. Ideally, the phase.dyna.tissueProfiles()
        #method would *ALWAYS* be implicitly called by the SimPhase.__init__()
        #method. Unfortunately, the non-trivial complexity of cell cluster
        #initialization requires we do so manually for now. Sad sandlion frowns!
        phase.dyna.tissueProfiles(sim, cells, self._p)

        # Display and/or save all simulation exports (e.g., animations).
        exppipe.pipeline(phase)

        #FIXME: Split each of the following blocks performing both plotting and
        #animating into their appropriate plotpipe.pipeline() or
        #animpipe.pipeline() functions.

        # run the molecules plots:
        if self._p.molecules_enabled and sim.molecules is not None:
            # reinit settings for plots, in case they've changed:
            sim.molecules.core.plot_init(self._p.network_config, self._p)
            sim.molecules.core.init_saving(cells, self._p, plot_type='sim')
            sim.molecules.core.export_all_data(sim, cells, self._p, message='auxiliary molecules')
            sim.molecules.core.plot(sim, cells, self._p, message='auxiliary molecules')
            sim.molecules.core.anim(phase=phase, message='auxiliary molecules')

        if self._p.grn_enabled and sim.grn is not None:
            # reinitialize the plot settings:
            sim.grn.core.plot_init(self._p.grn.conf, self._p)
            sim.grn.core.init_saving(cells, self._p, plot_type='sim', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, self._p, message='GRN molecules')
            sim.grn.core.plot(sim, cells, self._p, message='GRN molecules')
            sim.grn.core.anim(phase=phase, message='GRN molecules')

        #FIXME: Integrate into the plot pipeline.
        if self._p.Ca_dyn and self._p.ions_dict['Ca'] == 1:
            sim.endo_retic.init_saving(
                cells, self._p, plot_type='sim', nested_folder_name='ER')
            sim.endo_retic.plot_er(sim, cells, self._p)

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if self._p.plot.is_after_sim_show:
            plt.show()

        # Return this phase.
        return phase


    def plot_grn(self) -> SimPhase:
        '''
        Visualize the pure gene regulatory network (GRN) initialized and
        simulated by a prior call to the :meth:`sim_grn` method and export the
        resulting plots and animations to various output files, specified by the
        current configuration file.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # High-level simulation objects.
        sim = Simulator(self._p)

        #FIXME: See above for pertinent commentary. Mists of time, unpart!
        self._p.set_time_profile(phase_kind)  # force the time profile to be initialize

        # MoG = MasterOfGenes(self._p)
        MoG, cells, _ = fh.loadSim(self._p.grn_pickle_filename)

        # Simulation phase.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        sim.init_core(phase)

        # Initialize simulation data structures
        # sim.baseInit_all(cells, p)

        #FIXME: This... looks a bit hacky. That said, maybe it is good? Thunder
        #dragons unite!
        sim.time = MoG.time

        MoG.core.plot_init(self._p.grn.conf, self._p)
        MoG.core.init_saving(cells, self._p, plot_type='grn', nested_folder_name='RESULTS')
        MoG.core.export_all_data(sim, cells, self._p, message='gene products')
        MoG.core.plot(sim, cells, self._p, message='gene products')
        MoG.core.anim(phase=phase, message='gene products')

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if self._p.plot.is_after_sim_show:
            plt.show()

        # Return this phase.
        return phase

    # ..................{ EXCEPTIONS                         }..................
    @type_check
    def _die_if_seed_differs(
        self,
        p_old: Parameters,
        p_new: Parameters,
    ) -> None:
        '''
        Raise an exception if one or more general or seed (i.e., world) options
        differ between the passed simulation configurations, implying the
        current configuration to have been unsafely modified since the initial
        creation of the cell cluster for this configuration.

        Attributes
        ----------
        p_old : Parameters
            Previous simulation configuration to be compared.
        p_new : Parameters
            Current simulation configuration to be compared.

        Raises
        ----------
        :class:`BetseSimConfException`
            If these options differ for these configurations.
        '''

        if (p_old.config['general options'] != p_new.config['general options'] or
            p_old.config['world options'  ] != p_new.config['world options']):
            raise BetseSimConfException(
                'Important config file options are out of sync between '
                'seed and this init/sim attempt! '
                'Run "betse seed" again to match the current settings of '
                'this config file.')

    # ..................{ DEPRECATED                         }..................
    # The following methods have been deprecated for compliance with PEP 8.

    #FIXME: Remove all deprecated methods defined below *AFTER* a sufficient
    #amount of time -- say, mid to late 2018.

    @deprecated
    def makeWorld(self) -> None:
        return self.seed()

    @deprecated
    def initialize(self) -> None:
        return self.init()

    @deprecated
    def simulate(self) -> None:
        return self.sim()

    @deprecated
    def plotInit(self) -> None:
        return self.plot_init()

    @deprecated
    def plotSim(self) -> None:
        return self.plot_sim()
