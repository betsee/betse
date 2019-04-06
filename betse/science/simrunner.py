#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................
import matplotlib.pyplot as plt
import numpy as np
from betse.exceptions import BetseSimException, BetseSimConfException
from betse.lib.pickle import pickles
from betse.science import filehandling as fh
from betse.science.cells import Cells
from betse.science.chemistry.gene import MasterOfGenes
from betse.science.enum.enumconf import GrnUnpicklePhaseType
from betse.science.pipe.export.pipeexps import SimPipesExport
from betse.science.parameters import Parameters
from betse.science.phase import phasecallbacks
from betse.science.phase.phasecallbacks import SimCallbacksBCOrNoneTypes
from betse.science.phase.phasecls import SimPhase
from betse.science.enum.enumphase import SimPhaseKind
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.type.decorator.decorators import deprecated
from betse.util.type.decorator.decprof import log_time_seconds
from betse.util.type.text.string import strs
from betse.util.type.types import type_check
from matplotlib.collections import LineCollection, PolyCollection

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
    _callbacks : SimCallbacksBC
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
        callbacks: SimCallbacksBCOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this simulation runner.

        Attributes
        ----------
        p : Parameters
            Simulation configuration to be run. If this configuration is in the
            unloaded state (i.e., has *not* yet been loaded into memory from a
            YAML-formatted file on disk), an exception is raised.
        callbacks : SimCallbacksBCOrNoneTypes
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

    # ..................{ RUNNERS ~ seed                    }..................
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

        # Log this attempt.
        logs.log_info('Seeding simulation...')

        # Cumulative number of times that each call of this subcommand calls
        # the SimCallbacksBC.progressed() callback or a callback calling that
        # callback (e.g., SimCallbacksBC.progressed_next()).
        #
        # This magic number *must* be manually synchronized with the
        # implementation of both this subcommand and methods transitively
        # called by this subcommand (e.g., Cells.make_world()). Failure to do
        # so *will* result in fatal exceptions. Sadly, there exists no
        # reasonable means of either automating or enforcing this constraint.
        SEED_PROGRESS_TOTAL = (
            # Number of progress callbacks performed directly in this method.
            6 +
            # Number of progress callbacks performed by Cells.make_world().
            Cells.MAKE_WORLD_PROGRESS_TOTAL
        )

        # Notify the caller of the range of work performed by this subcommand.
        self._callbacks.progress_ranged(progress_max=SEED_PROGRESS_TOTAL)

        # Simulation phase.
        phase = SimPhase(
            kind=SimPhaseKind.SEED, p=self._p, callbacks=self._callbacks)

        # Create the pseudo-randomized cell cluster.
        phase.cells.make_world(phase)

        # Initialize core simulation data structures.
        self._callbacks.progressed_next(
            status='Creating core computational matrices...')
        phase.sim.init_core(phase)

        # Define the tissue and boundary profiles for plotting.
        self._callbacks.progressed_next(
            status='Creating tissue, surgery, and boundary profiles...')
        phase.dyna.init_profiles(phase)

        # Redo gap junctions to isolate different tissue types.
        self._callbacks.progressed_next(
            status='Creating gap junction connection network...')
        phase.cells.redo_gj(phase)

        # Create a Laplacian and solver for discrete transfers on closed,
        # irregular cell network.
        self._callbacks.progressed_next(
            status='Creating cell network Poisson solver...')
        phase.cells.graphLaplacian(self._p)

        # Create accessory matrices depending on user requirements.
        if self._p.deformation:
            phase.cells.deform_tools(self._p)

        # if self._p.sim_eosmosis:
        #     phase.cells.eosmo_tools(self._p)

        # Pickle this cell cluster to disk.
        self._callbacks.progressed_next(
            status='Saving seeded cell cluster...')
        phase.cells.save_cluster(phase)

        # Log the completion of this phase.
        logs.log_info('Cell cluster creation complete!')
        phase.sim.sim_info_report(phase)

        # Signal the completion of this phase with respect to progress.
        self._callbacks.progressed_last()

        # Return this phase.
        return phase

    # ..................{ RUNNERS ~ (init|sim)              }..................
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

        if not files.is_file(self._p.seed_pickle_filename):
            if not self._p.autoInit:
                #FIXME: Call the _die_unless_file_pickled() function instead
                #both here and everywhere else we perform similar logic. Bugbears!
                raise BetseSimException(
                    'Initialization halted due to missing seed. '
                    'Please run "betse seed" and try again.')

            # Create an instance of world.
            logs.log_info('Automatically seeding cell cluster...')
            self.seed()
            logs.log_info('Now using seed to run initialization.')

        # Load the seed from cache.
        cells, p_old = fh.loadWorld(self._p.seed_pickle_filename)
        logs.log_info('Cell cluster loaded.')

        # check to ensure compatibility between original and present sim files:
        self._die_if_seed_differs(p_old, self._p)

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        phase.sim.init_core(phase)

        # Run this simulation phase.
        phase.sim.sim_info_report(phase)
        phase.sim.run_sim_core(phase)

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

        if not files.is_file(self._p.init_pickle_filename):
            if not self._p.autoInit:
                raise BetseSimException(
                    'Simulation halted due to missing initialization. '
                    'Please run "betse init" and try again.')

            logs.log_info('Automatically initializing cell cluster...')
            self.init()
            logs.log_info('Now using initialization to run simulation.')

        # Load the initialization from cache.
        sim, cells, p_old = fh.loadSim(self._p.init_pickle_filename)

        # Ensure compatibility between original and present config files.
        self._die_if_seed_differs(p_old, self._p)

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        # Run and save the simulation to the cache.
        sim.sim_info_report(phase)
        sim.run_sim_core(phase)

        # Return this phase.
        return phase

    # ..................{ RUNNERS ~ grn                     }..................
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

        # Simulator objects initialized below.
        cells = None
        sim = None

        # Log this simulation.
        logs.log_info(
            'Running gene regulatory network "%s" '
            'defined in config file "%s"...',
            pathnames.get_basename(self._p.grn_config_filename),
            self._p.conf_basename)

        # If networking an uninitialized, unsimulated cell cluster...
        if self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.SEED:
            if not files.is_file(self._p.seed_pickle_filename):
                if not self._p.autoInit:
                    raise BetseSimException(
                        'Simulation halted due to missing core seed. '
                        'Please run "betse seed" and try again.')

                # Create an instance of world.
                logs.log_info('Automatically seeding cell cluster...')
                self.seed()

            # Load the seed from cache.
            cells, _ = fh.loadWorld(self._p.seed_pickle_filename)
            logs.log_info('Running gene regulatory network on betse seed...')
            logs.log_info('Now using cell cluster to run initialization.')

            # Simulation phase.
            phase = SimPhase(
                kind=phase_kind,
                callbacks=self._callbacks,
                cells=cells,
                p=self._p,
            )

            # Initialize core simulation data structures.
            phase.sim.init_core(phase)
            phase.sim.init_dynamics(phase)

            #FIXME: Shift the following assignments into a new public
            #"Simulator" method -- say, Simulator.init_core_null().

            # Initialize other aspects required for piggyback of GRN on the
            # sim object.
            phase.sim.time = []
            phase.sim.vm = -50e-3 * np.ones(phase.sim.mdl)

            # Initialize key fields of simulator required to interface
            # (dummy init).
            phase.sim.rho_pump = 1.0
            phase.sim.rho_channel = 1.0
            phase.sim.conc_J_x = np.zeros(phase.sim.edl)
            phase.sim.conc_J_y = np.zeros(phase.sim.edl)
            phase.sim.J_env_x  = np.zeros(phase.sim.edl)
            phase.sim.J_env_y  = np.zeros(phase.sim.edl)
            phase.sim.u_env_x  = np.zeros(phase.sim.edl)
            phase.sim.u_env_y  = np.zeros(phase.sim.edl)
        # Else if networking an initialized but unsimulated cell cluster...
        elif self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.INIT:
            if not files.is_file(self._p.init_pickle_filename):
                if not self._p.autoInit:
                    raise BetseSimException(
                        'Simulation halted due to missing core initialization. '
                        'Please run "betse init" and try again.')

                logs.log_info('Automatically initializing cell cluster...')
                self.init()
                logs.log_info('Now using initialization to run simulation.')

            # Load the initialization from cache.
            logs.log_info('Running gene regulatory network on betse init...')
            sim, cells, _ = fh.loadSim(self._p.init_pickle_filename)
        # Else if networking an initialized, simulated cell cluster...
        elif self._p.grn_unpickle_phase_type is GrnUnpicklePhaseType.SIM:
            if not files.is_file(self._p.sim_pickle_filename):
                raise BetseSimException(
                    'Simulation halted due to missing core simulation. '
                    'Please run "betse sim" and try again.')

            # Load the simulation from cache.
            logs.log_info('Running gene regulatory network on betse sim...')
            sim, cells, _ = fh.loadSim(self._p.sim_pickle_filename)
        # Else, this type of networking is unrecognized. Raise an exception.
        else:
            raise BetseSimConfException(
                'Gene regulatory network (GRN) unpickle simulation phase '
                '"{}" unrecognized.'.format(self._p.grn_unpickle_phase_type))

        # If *NOT* defined above, define this simulation phase.
        if phase is None:
            phase = SimPhase(
                kind=phase_kind,
                callbacks=self._callbacks,
                cells=cells,
                p=self._p,
                sim=sim,
            )

            # Reinitialize all profiles.
            phase.dyna.init_profiles(phase)
            phase.dyna.init_events(phase)

        # If *NOT* restarting from a prior GRN run, start a new GRN.
        if self._p.grn_unpickle_filename is None:
            # Log this start.
            logs.log_info("Initializing the gene regulatory network...")

            # Create and initialize an instance of master of metabolism.
            MoG = MasterOfGenes(self._p)
            MoG.read_gene_config(phase)
        # Else, restart from a prior GRN run.
        else:
            # Log this restart.
            logs.log_info(
                'Reinitializing the gene regulatory network from "%s"...',
                pathnames.get_basename(self._p.grn_unpickle_filename))

            # If this file does *NOT* exist, raise an exception.
            _die_unless_file_pickled(
                filename=self._p.grn_unpickle_filename,
                subcommand='sim-grn',
                subcommand_label='Gene regulatory network')

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

                # If no prior initialization exists, raise an exception.
                if not files.is_file(self._p.sim_pickle_filename):
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
                sim_old, cells_old, _ = fh.loadSim(self._p.sim_pickle_filename)

                #FIXME: This phase object would ideally be pickled to and
                #from the "self._p.sim_pickle_filename" file loaded above, in
                #which case this local variable would be safely removable.

                # Original simulation phase. To avoid caller confusion, the
                # optional "callbacks" parameter is intentionally *NOT* passed.
                phase_old = SimPhase(
                    kind=phase_kind,
                    p=self._p,
                    cells=cells_old,
                    sim=sim_old,
                )

                # Initialize all tissue profiles on original cells.
                phase_old.dyna.init_profiles(phase_old)

                for cut_profile_name in phase_old.p.event_cut_profile_names:
                    logs.log_info(
                        'Cutting cell cluster via cut profile "%s"...',
                        cut_profile_name)

                    # Object picking the cells removed by this cut profile.
                    tissue_picker = phase_old.dyna.cut_name_to_profile[
                        cut_profile_name].picker

                    # One-dimensional Numpy arrays of the indices of all
                    # cells and cell membranes to be removed.
                    target_inds_cell, target_inds_mems = (
                        tissue_picker.pick_cells_and_mems(
                            cells=cells_old, p=self._p))

                    MoG.core.mod_after_cut_event(
                        phase, target_inds_cell, target_inds_mems)
                    MoG.core.redefine_dynamic_dics(sim, cells, self._p)

                    logs.log_info(
                        'Reinitializing gene regulatory network '
                        'for simulation...')
                    MoG.reinitialize(phase)

            # if self._p.use_microtubules:
            #     sim.mtubes.mtubes_x = MoG.mtubes_x_time[-1]
            #     sim.mtubes.mtubes_y = MoG.mtubes_y_time[-1]
            #
            #     sim.mtubes.uxmt, sim.mtubes.uymt = sim.mtubes.mtubes_to_cell(
            #         cells, self._p)

        logs.log_info('Simulating gene regulatory network...')
        MoG.run_core_sim(phase)

        # Return this phase.
        return phase

    # ..................{ PLOTTERS                          }..................
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

        # Log this plotting attempt.
        logs.log_info(
            'Plotting cell cluster with configuration file "%s".',
            self._p.conf_basename)

        # If an initialization does *NOT* already exist, raise an exception.
        _die_unless_file_pickled(
            filename=self._p.seed_pickle_filename,
            subcommand='seed',
            subcommand_label='Seed')

        # Load the seed from cache.
        cells, _ = fh.loadWorld(self._p.seed_pickle_filename)
        logs.log_info('Cell cluster loaded.')

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=SimPhaseKind.SEED,
            p=self._p,
            cells=cells,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        phase.sim.init_core(phase)
        phase.dyna.init_profiles(phase)

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
                np.log10(phase.sim.D_env_weight.reshape(phase.cells.X.shape)),
                origin='lower',
                extent=[
                    self._p.um * phase.cells.xmin,
                    self._p.um * phase.cells.xmax,
                    self._p.um * phase.cells.ymin,
                    self._p.um * phase.cells.ymax,
                ],
                cmap=self._p.background_cm,
            )
            plt.colorbar()

            cell_edges_flat = self._p.um * phase.cells.mem_edges_flat
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
                extent=[
                    self._p.um * phase.cells.xmin,
                    self._p.um * phase.cells.xmax,
                    self._p.um * phase.cells.ymin,
                    self._p.um * phase.cells.ymax,
                ],
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
            base_points = np.multiply(phase.cells.cell_verts, self._p.um)
            col_cells = PolyCollection(
                base_points, facecolors='k', edgecolors='none')
            col_cells.set_alpha(0.3)
            ax_x.add_collection(col_cells)

        con_segs = phase.cells.nn_edges
        connects = self._p.um * np.asarray(con_segs)
        collection = LineCollection(connects, linewidths=1.0, color='b')
        ax_x.add_collection(collection)
        plt.axis('equal')
        plt.axis([
            self._p.um * phase.cells.xmin,
            self._p.um * phase.cells.xmax,
            self._p.um * phase.cells.ymin,
            self._p.um * phase.cells.ymax,
        ])

        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Cell Connectivity Network')

        if self._p.autosave:
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

        # If an initialization does *NOT* already exist, raise an exception.
        _die_unless_file_pickled(
            filename=self._p.init_pickle_filename,
            subcommand='init',
            subcommand_label='Initialization')

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # Load the initialization from cache.
        sim, cells, _ = fh.loadSim(self._p.init_pickle_filename)

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        #FIXME: This... isn't the best. Ideally, the phase.dyna.init_profiles()
        #method would *ALWAYS* be implicitly called by the SimPhase.__init__()
        #method. Unfortunately, the non-trivial complexity of cell cluster
        #initialization requires we do so manually for now. Sad sandlion frowns!
        phase.dyna.init_profiles(phase)

        # Display and/or save all initialization exports (e.g., animations).
        SimPipesExport().export(phase)

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

        # If a simulation does *NOT* already exist, raise an exception.
        _die_unless_file_pickled(
            filename=self._p.sim_pickle_filename,
            subcommand='sim',
            subcommand_label='Simulation')

        # Simulation phase type.
        phase_kind = SimPhaseKind.SIM

        # Load the simulation from the cache.
        sim, cells, _ = fh.loadSim(self._p.sim_pickle_filename)

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            sim=sim,
            callbacks=self._callbacks,
        )

        #FIXME: This... isn't the best. Ideally, the phase.dyna.init_profiles()
        #method would *ALWAYS* be implicitly called by the SimPhase.__init__()
        #method. Unfortunately, the non-trivial complexity of cell cluster
        #initialization requires we do so manually for now. Sad sandlion frowns!
        phase.dyna.init_profiles(phase)

        # Display and/or save all simulation exports (e.g., animations).
        SimPipesExport().export(phase)

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

        # MoG = MasterOfGenes(self._p)
        MoG, cells, _ = fh.loadSim(self._p.grn_pickle_filename)

        # Simulation phase.
        phase = SimPhase(
            kind=phase_kind,
            cells=cells,
            p=self._p,
            callbacks=self._callbacks,
        )

        # Initialize core simulation data structures.
        phase.sim.init_core(phase)

        # Initialize simulation data structures
        # sim.baseInit_all(cells, p)

        #FIXME: This... looks a bit hacky. That said, maybe it is good? Thunder
        #dragons unite!
        phase.sim.time = MoG.time

        MoG.core.plot_init(self._p.grn.conf, self._p)
        MoG.core.init_saving(
            cells, self._p, plot_type='grn', nested_folder_name='RESULTS')
        MoG.core.export_all_data(
            phase.sim, cells, self._p, message='gene products')
        MoG.core.plot(
            phase.sim, cells, self._p, message='gene products')
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

    # ..................{ DEPRECATED                        }..................
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

# ....................{ EXCEPTIONS                        }....................
@type_check
def _die_unless_file_pickled(
    filename: str, subcommand: str, subcommand_label: str) -> None:
    '''
    Raise an exception containing the passed machine-readable subcommand and
    human-readable label denoting that subcommand unless the file with the
    passed filename pickled by a prior run of that subcommand exists.

    Parameters
    ----------
    filename : str
        Absolute filename of the pickled file to be validated.
    subcommand : str
        Machine-readable case-sensitive name of the prior subcommand expected
        to have pickled this file (e.g., ``seed``, ``init``).
    subcommand_label : str
        Human-readable case-insensitive noun denoting the type of prior
        subcommand expected to have pickled this file (e.g., ``Seed``,
        ``initialization``).

    Raises
    ----------
    BetseSimException
        If this file does *not* exist.
    '''

    # If this file does *NOT* exist...
    if not files.is_file(filename):
        # Uppercase the first character of this label for readability.
        subcommand_label_cased = strs.uppercase_char_first(subcommand_label)

        # Raise an exception embedding these parameters.
        raise BetseSimException(
            '{} not previously run; '
            'please run "betse {}" and try again.\n'
            'Specifically, file not found: {}'.format(
                subcommand_label_cased, subcommand, filename))
