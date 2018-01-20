#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import time
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from betse.exceptions import (
    BetseFileException, BetseSimException, BetseSimConfigException)
from betse.science import filehandling as fh
from betse.science.cells import Cells
from betse.science.chemistry.gene import MasterOfGenes
# from betse.science.math import toolbox as tb
from betse.science.config import confio
from betse.science.export import exppipe
from betse.science.parameters import Parameters
from betse.science.sim import Simulator
from betse.science.simulate.simphase import SimPhase, SimPhaseKind
from betse.science.tissue.tishandler import TissueHandler
from betse.science.visual.plot import plotutil as viz
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.type.call.callables import deprecated
from betse.lib.pickle import pickles

# ....................{ CLASSES                            }....................
class SimRunner(object):
    '''
    High-level simulation class encapsulating the running of *all* available
    simulation phases.

    This class provides high-level methods for initializing, running, and
    plotting simulations specified by the YAML configuration file with which
    this class is instantiated. Thus, each instance of this class only handles a
    single simulation.

    Attributes
    ----------
    _config_filename : str
        Absolute path of the YAML file configuring this simulation.
    _config_basename : str
        Basename of the YAML file configuring this simulation.
    '''

    def __init__(self, conf_filename: str) -> None:

        super().__init__()

        # Validate and localize this filename.
        files.die_unless_file(conf_filename)
        self._config_filename = conf_filename
        self._config_basename = pathnames.get_basename(self._config_filename)

    # ..................{ RUNNERS                            }..................
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

        # Simulation configuration, simulator, and cell cluster.
        p = Parameters().load(self._config_filename)
        sim = Simulator(p)
        cells = Cells(p)

        # Simulation phase.
        phase = SimPhase(
            kind=SimPhaseKind.SEED, cells=cells, p=p, sim=sim)

        logs.log_info('Creating cell cluster...')
        cells.make_world(phase)  # call function to create the world

        # define the tissue and boundary profiles for plotting:
        logs.log_info('Defining tissue and boundary profiles...')
        sim.baseInit_all(cells, p)
        dyna = TissueHandler(sim, cells, p)
        dyna.tissueProfiles(sim, cells, p)
        cells.redo_gj(dyna, p)  # redo gap junctions to isolate different tissue types

        # make a laplacian and solver for discrete transfers on closed, irregular cell network
        logs.log_info('Creating cell network Poisson solver...')
        cells.graphLaplacian(p)

        if p.td_deform is False:  # if time-dependent deformation is not required
            cells.lapGJ = None
            cells.lapGJ_P = None  # null out the non-inverse matrices -- we don't need them

        # make accessory matrices depending on user requirements:
        if p.fluid_flow is True or p.deformation is True:
            if p.deformation is True:
                cells.deform_tools(p)

        # if p.sim_eosmosis is True:
        #     cells.eosmo_tools(p)

        # finish up:
        cells.save_cluster(p)
        logs.log_info('Cell cluster creation complete!')

        sim.sim_info_report(cells,p)

        # Return this phase.
        return phase

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

        start_time = time.time()  # get a start value for timing the simulation

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        #FIXME: The Parameters.__init__() method should *REQUIRE* that a time
        #profile type be passed. The current approach leaves critical attributes
        #undefined in the event that the optional Parameters.set_time_profile()
        #method is left uncalled, which is pretty unacceptable.
        #FIXME: Actually, no. All logic performed by the set_time_profile()
        #method should be shifted into the SimPhase.__init__() method. See a
        #FIXME comment preceding the set_time_profile() method for details.

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be initialize
        p.run_sim = False # let the simulator know we're just running an initialization

        # Simulation cell cluster.
        # cells, _ = fh.loadSim(cells.savedWorld)
        cells = Cells(p)  # create an instance of world

        #FIXME: This if conditional is repeated verbatim twice below. Generalize
        #into a new _load_cells() method containing this if conditional and
        #returning the loaded "Cells" instance; then, call this method both here
        #and everywhere this repeated logic appears below. Starbust dragons!

        if files.is_file(cells.savedWorld):
            cells,p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')

            # check to ensure compatibility between original and present sim files:
            self._die_unless_seed_same(p_old, p)

        else:
            logs.log_warning("Ooops! No such cell cluster file found to load!")

            if p.autoInit:
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

        # Simulation simulator.
        sim = Simulator(p=p)

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

        sim.run_sim = False

        # Initialize simulation data structures, run, and save simulation phase
        sim.baseInit_all(cells, p)
        sim.sim_info_report(cells, p)
        sim.run_sim_core(phase)

        logs.log_info(
            'Initialization completed in %d seconds.',
            round(time.time() - start_time, 2))

        # Return this phase.
        return phase

    def sim(self) -> SimPhase:
        '''
        Simulate this simulation with the cell cluster initialized by a prior
        call to the :meth:`init` method and cache this simulation to an output
        file, specified by the current configuration file.

        This method _must_ be called prior to the :meth:`:meth:`plot_sim`
        method, which consumes this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        # Log this attempt.
        logs.log_info('Running simulation...')

        start_time = time.time()  # get a start value for timing the simulation

        # Simulation phase type.
        phase_kind = SimPhaseKind.SIM

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be initialize
        p.run_sim = True    # set on the fly a boolean to let simulator know we're running a full simulation

        # Simulation simulator.
        sim = Simulator(p=p)

        if files.is_file(sim.savedInit):
            sim,cells, p_old = fh.loadSim(sim.savedInit)  # load the initialization from cache

            # check to ensure compatibility between original and present sim files:
            self._die_unless_seed_same(p_old, p)

        else:
            logs.log_warning(
                "No initialization file found to run this simulation!")

            if p.autoInit:
                logs.log_info("Automatically running initialization...")
                self.init()
                logs.log_info('Now using initialization to run simulation.')
                sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

            else:
                raise BetseSimException(
                    'Simulation terminated due to missing initialization. '
                    'Please run an initialization and try again.')

        # Simulation phase, created *AFTER* unpickling these objects above.
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

        # Reinitialize save and load directories in case params defines new ones
        # for this sim.
        sim.fileInit(p)

        # Run and save the simulation to the cache.
        sim.sim_info_report(cells, p)
        sim.run_sim_core(phase)

        logs.log_info(
            'Simulation completed in %d seconds.',
            round(time.time() - start_time, 2))

        # Return this phase.
        return phase


    def sim_grn(self) -> SimPhase:
        '''
        Initialize and simulate a pure gene regulatory network (GRN) _without_
        bioelectrics with the cell cluster seeded by a prior call to the
        :meth:`seed` method and cache this initialization and simulation to
        output files, specified by the current configuration file.

        This method _must_ be called prior to the :meth:`plot_grn` method, which
        consumes this output as input.

        Returns
        ----------
        SimPhase
            High-level simulation phase instance encapsulating all objects
            internally created by this method to run this phase.
        '''

        start_time = time.time()  # get a start value for timing the simulation

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be initialize
        p.run_sim = False

        logs.log_info(
            ('Running gene regulatory network {} defined in config file {}.').format(p.grn_config_filename, self._config_basename))

        # cells object:
        cells = Cells(p)

        # Simulation simulator.
        sim = Simulator(p=p)

        if p.grn_piggyback == 'seed':

            if files.is_file(cells.savedWorld):
                cells, p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
                logs.log_info('Running gene regulatory network on betse seed...')

                # Initialize simulation data structures
                sim.baseInit_all(cells, p)

                # Initialize other aspects required for piggyback of GRN on the sim object:
                sim.time = []
                sim.vm = -50e-3 * np.ones(sim.mdl)
                # initialize key fields of simulator required to interface (dummy init)
                sim.rho_pump = 1.0
                sim.rho_channel = 1.0
                sim.conc_J_x = np.zeros(sim.edl)
                sim.conc_J_y = np.zeros(sim.edl)
                sim.J_env_x = np.zeros(sim.edl)
                sim.J_env_y = np.zeros(sim.edl)
                sim.u_env_x = np.zeros(sim.edl)
                sim.u_env_y = np.zeros(sim.edl)

            else:
                logs.log_warning("Ooops! No such cell cluster file found to load!")

                if p.autoInit:
                    logs.log_info(
                        'Automatically seeding cell cluster from config file settings...')
                    self.seed()  # create an instance of world
                    logs.log_info(
                        'Now using cell cluster to run initialization.')
                    cells, _ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

                else:
                    raise BetseSimException(
                        "Run terminated due to missing seed.\n"
                        "Please run 'betse seed' to try again.")

        elif p.grn_piggyback == 'init':
            if files.is_file(sim.savedInit):
                logs.log_info('Running gene regulatory network on betse init...')
                sim, cells, p_old = fh.loadSim(sim.savedInit)  # load the initialization from cache

            else:
                logs.log_warning(
                    "No initialization file found to run the GRN simulation!")

                if p.autoInit:
                    logs.log_info("Automatically running initialization...")
                    self.init()
                    logs.log_info('Now using initialization to run simulation.')
                    sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

                else:
                    raise BetseSimException(
                        'Simulation terminated due to missing core initialization. '
                        'Please run a betse initialization and try again.')

        elif p.grn_piggyback == 'sim':
            if files.is_file(sim.savedSim):
                logs.log_info('Running gene regulatory network on betse sim...')
                sim, cells, p_old = fh.loadSim(sim.savedSim)  # load the initialization from cache

            else:
                logs.log_warning(
                    "No simulation file found to run the GRN simulation!")
                raise BetseSimException(
                    'Simulation terminated due to missing core simulation. '
                    'Please run a betse simulation and try again.')

        # # Simulation simulator.
        # sim = Simulator(p=p)

        # Simulation phase.
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

        # If loading from a previously pickled "sim-grn" file...
        if p.loadMoG is not None and files.is_file(p.loadMoG):
            # load previously run instance of master of genes:
            MoG, _, _ = pickles.load(p.loadMoG)

            logs.log_info(("Reinitializing the gene regulatory network from {}...").format(p.loadMoG))

            is_cut_done = sim.dyna.event_cut.is_fired

            # if running on a sim with a cut event, must remove cells:
            if sim.dyna.event_cut is not None and is_cut_done:
                simu = Simulator(p=p)

                logs.log_info(
                    'A cutting event has been run, so the GRN object needs to be modified...')

                if files.is_file(simu.savedInit):
                    logs.log_info(
                        'Loading betse init from cache for reference to original cells...')

                    init, cellso, p_old = fh.loadSim(simu.savedInit)  # load the initialization from cache

                    dyna = TissueHandler(init, cellso, p)  # create the tissue dynamics object on original cells
                    dyna.tissueProfiles(init, cellso, p)  # initialize all tissue profiles on original cells

                    for cut_profile_name in dyna.event_cut.profile_names:
                        logs.log_info(
                            'Cutting cell cluster via cut profile "%s"...',
                            cut_profile_name)

                        # Object picking the cells removed by this cut profile.
                        tissue_picker = dyna.cut_name_to_profile[
                            cut_profile_name].picker

                        # One-dimensional Numpy arrays of the indices of all
                        # cells and cell membranes to be removed.
                        target_inds_cell, target_inds_mem = (
                            tissue_picker.pick_cells_and_mems(
                                cells=cellso, p=p))

                        MoG.core.mod_after_cut_event(
                            target_inds_cell, target_inds_mem, sim, cells, p)

                        logs.log_info(
                            "Redefining dynamic dictionaries to point to the new sim...")
                        MoG.core.redefine_dynamic_dics(sim, cells, p)

                        logs.log_info(
                            "Reinitializing the gene regulatory network for simulation...")
                        MoG.reinitialize(sim, cells, p)

                else:
                    logs.log_warning(
                        "This situation is complex due to a cutting event being run.\n"
                        "Please have a corresponding init file to run the GRN simulation!")

                    raise BetseSimException(
                        'Simulation terminated due to missing core init. '
                        'Please alter GRN settings and try again.')

        # Else, a previously pickled "sim-grn" file is *NOT* being loaded from.
        # In this case, create a new GRN from scratch.
        else:
            # If GRN support is disabled, raise an exception.
            if not p.grn_enabled:
                raise BetseSimException('GRN support disabled.')

            # create an instance of master of metabolism
            MoG = MasterOfGenes(p)

            # initialize it:
            logs.log_info("Initializing the gene regulatory network...")
            MoG.read_gene_config(sim, cells, p)

        logs.log_info("Running gene regulatory network test simulation...")

        MoG.run_core_sim(sim, cells, p)

        logs.log_info(
            'Gene regulatory network test completed in %d seconds.',
            round(time.time() - start_time, 2))

        # Return this phase.
        return phase

    # ..................{ PLOTTERS                           }..................
    #FIXME: Shift the low-level matplotlib plotting performed by this method
    #into a new "betse.science.visual.seedpipe" submodule.

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
            self._config_basename)

        # Simulation configuration, simulator, and cell cluster.
        p = Parameters().load(self._config_filename)
        sim = Simulator(p)
        cells = Cells(p)

        if files.is_file(cells.savedWorld):
            cells, _ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')
        else:
            raise BetseSimException(
                "Ooops! No such cell cluster file found to load!")

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(
            kind=SimPhaseKind.SEED, cells=cells, p=p, sim=sim)

        sim.baseInit_all(cells,p)
        dyna = TissueHandler(sim,cells,p)
        dyna.tissueProfiles(sim,cells,p)

        #FIXME: Shift everything below into a new seed-specific pipeline -- say,
        #betse.science.export.exppipe.pipeline_seed().
        if p.autosave:
            savedImg = pathnames.join(p.init_export_dirname, 'fig_')

        if p.plot_cell_cluster:
            fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(
                p, dyna, cells, clrmap=p.background_cm)

            if p.autosave:
                savename10 = savedImg + 'cluster_mosaic' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.plot.is_after_sim_show:
                plt.show(block = False)

        if p.is_ecm and p.plot_cluster_mask:
            plt.figure()
            ax99 = plt.subplot(111)
            plt.imshow(
                np.log10(sim.D_env_weight.reshape(cells.X.shape)),
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.background_cm,
            )
            plt.colorbar()

            cell_edges_flat = p.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(1.0)
            ax99.add_collection(coll)

            plt.title('Logarithm of Environmental Diffusion Weight Matrix')

            if p.autosave:
                savename10 = savedImg + 'env_diffusion_weights' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.plot.is_after_sim_show:
                plt.show(block = False)

            plt.figure()
            plt.imshow(
                cells.maskM,
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.background_cm,
            )
            plt.colorbar()
            plt.title('Cluster Masking Matrix')

            if p.autosave:
                savename = savedImg + 'cluster_mask' + '.png'
                plt.savefig(savename,format='png',transparent=True)

            if p.plot.is_after_sim_show:
                plt.show(block = False)

        # Plot gap junctions.
        if p.plot_cell_connectivity:
            plt.figure()
            ax_x = plt.subplot(111)

            if p.showCells:
                base_points = np.multiply(cells.cell_verts, p.um)
                col_cells = PolyCollection(base_points, facecolors='k', edgecolors='none')
                col_cells.set_alpha(0.3)
                ax_x.add_collection(col_cells)

            con_segs = cells.nn_edges
            connects = p.um*np.asarray(con_segs)
            collection = LineCollection(connects,linewidths=1.0,color='b')
            ax_x.add_collection(collection)
            plt.axis('equal')
            plt.axis([cells.xmin*p.um,cells.xmax*p.um,cells.ymin*p.um,cells.ymax*p.um])

            ax_x.set_xlabel('Spatial x [um]')
            ax_x.set_ylabel('Spatial y [um')
            ax_x.set_title('Cell Connectivity Network')

            if p.autosave is True:
                savename10 = savedImg + 'gj_connectivity_network' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        if p.turn_all_plots_off is False:
            plt.show(block=False)
            plt.show()

        else:
            logs.log_info(
                'Plots exported to init results folder '
                'defined in configuration file "%s".',
                self._config_basename)

        # Return this phase.
        return phase

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
            self._config_basename)

        # Simulation phase type.
        phase_kind = SimPhaseKind.INIT

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be initialize

        # Simulation simulator.
        sim = Simulator(p=p)

        #FIXME: Bizarre logic. We create a "Simulator" instance above only to
        #test whether a single file exists and, if so, replace that instance
        #with a pickled "Simulator" instance unpickled from that file. Let's cut
        #out the API middleman, as it were, by obtaining the pathname for this
        #file from the "Parameters" instance instead and then removing the above
        #instantiation of "sim = Simulator(p)".

        if files.is_file(sim.savedInit):
            sim, cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        else:
            raise BetseSimException(
                "Ooops! No such initialization file found to plot!")

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

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
        if p.molecules_enabled and sim.molecules is not None:
            # reinit settings for plots, in case they've changed:
            sim.molecules.core.plot_init(p.network_config, p)
            sim.molecules.core.init_saving(cells, p, plot_type='init')
            sim.molecules.core.export_all_data(sim, cells, p, message='auxiliary molecules')
            sim.molecules.core.plot(sim, cells, p, message='auxiliary molecules')
            sim.molecules.core.anim(phase=phase, message='auxiliary molecules')

        if p.grn_enabled and sim.grn is not None:
            # reinitialize the plot settings:
            sim.grn.core.plot_init(p.grn.conf, p)
            sim.grn.core.init_saving(
                cells, p, plot_type='init', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, p, message='GRN molecules')
            sim.grn.core.plot(sim, cells, p, message='GRN molecules')
            sim.grn.core.anim(phase=phase, message='GRN molecules')

        #FIXME: Integrate into the plot pipeline.
        if p.Ca_dyn and p.ions_dict['Ca'] == 1:
            sim.endo_retic.init_saving(
                cells, p, plot_type='init', nested_folder_name='ER')
            sim.endo_retic.plot_er(sim, cells, p)

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if p.plot.is_after_sim_show:
            plt.show()

        # Return this phase.
        return phase

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
            self._config_basename)

        # Simulation phase type.
        phase_kind = SimPhaseKind.SIM

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be simulation

        # Simulation simulator.
        sim = Simulator(p=p)

        # If this simulation has yet to be run, fail.
        if not files.is_file(sim.savedSim):
            raise BetseFileException(
                'Simulation cache file "{}" not found to plot '
                '(e.g., due to no simulation having been run).'.format(
                    sim.savedSim))

        # Load the simulation from the cache.
        sim, cells, _ = fh.loadSim(sim.savedSim)

        # Simulation phase, created *AFTER* unpickling these objects above
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

        # Display and/or save all simulation exports (e.g., animations).
        exppipe.pipeline(phase)

        #FIXME: Split each of the following blocks performing both plotting and
        #animating into their appropriate plotpipe.pipeline() or
        #animpipe.pipeline() functions.

        # run the molecules plots:
        if p.molecules_enabled and sim.molecules is not None:
            # reinit settings for plots, in case they've changed:
            sim.molecules.core.plot_init(p.network_config, p)
            sim.molecules.core.init_saving(cells, p, plot_type='sim')
            sim.molecules.core.export_all_data(sim, cells, p, message='auxiliary molecules')
            sim.molecules.core.plot(sim, cells, p, message='auxiliary molecules')
            sim.molecules.core.anim(phase=phase, message='auxiliary molecules')

        if p.grn_enabled and sim.grn is not None:
            # reinitialize the plot settings:
            sim.grn.core.plot_init(p.grn.conf, p)
            sim.grn.core.init_saving(cells, p, plot_type='sim', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, p, message='GRN molecules')
            sim.grn.core.plot(sim, cells, p, message='GRN molecules')
            sim.grn.core.anim(phase=phase, message='GRN molecules')

        #FIXME: Integrate into the plot pipeline.
        if p.Ca_dyn and p.ions_dict['Ca'] == 1:
            sim.endo_retic.init_saving(
                cells, p, plot_type='sim', nested_folder_name='ER')
            sim.endo_retic.plot_er(sim, cells, p)

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if p.plot.is_after_sim_show:
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

        # Simulation configuration.
        p = Parameters().load(self._config_filename)
        p.set_time_profile(phase_kind)  # force the time profile to be initialize

        # MoG = MasterOfGenes(p)
        MoG, cells, _ = fh.loadSim(p.savedMoG)

        # Simulation simulator.
        sim = Simulator(p)

        # Simulation phase.
        phase = SimPhase(kind=phase_kind, cells=cells, p=p, sim=sim)

        # Initialize simulation data structures
        sim.baseInit_all(cells, p)
        sim.time = MoG.time

        MoG.core.plot_init(p.grn.conf, p)
        MoG.core.init_saving(cells, p, plot_type='grn', nested_folder_name='RESULTS')
        MoG.core.export_all_data(sim, cells, p, message='gene products')
        MoG.core.plot(sim, cells, p, message='gene products')
        MoG.core.anim(phase=phase, message='gene products')

        # If displaying plots, block on all previously plots previously
        # displayed as non-blocking. If this is *NOT* done, these plots will
        # effectively *NEVER* displayed be on most systems.
        if p.plot.is_after_sim_show:
            plt.show()

        # Return this phase.
        return phase

    # ..................{ UTILITIES                          }..................
    def _die_unless_seed_same(self, p_old, p) -> None:
        '''
        Raise an exception unless the two passed simulation configurations share
        the same general and seed (i.e., world) options, implying the current
        configuration to have been modified since the initial seeding of this
        configuration's cell cluster.
        '''

        if (p_old.config['general options'] != p.config['general options'] or
            p_old.config['world options'  ] != p.config['world options']):
            # logs.log_warning('---------------------------------------------------')
            # logs.log_warning('**WARNING!**')
            # logs.log_warning('Important config file options are out of sync ')
            # logs.log_warning('between the seed and this init/sim attempt! ')
            # logs.log_warning('Run "betse seed" again to match the current settings')
            # logs.log_warning(' of this config file.')
            # logs.log_warning('---------------------------------------------------')

            raise BetseSimConfigException(
                'Important config file options are out of sync between '
                'seed and this init/sim attempt! '
                'Run "betse seed" again to match the current settings of '
                'this config file.')

    # ..................{ DEPRECATED                         }..................
    # The following methods have been deprecated for compliance with PEP 8.

    #FIXME: Remove all deprecated methods defined below *AFTER* a sufficient
    #amount of time -- say, mid to late 2017.

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
    def plotWorld(self) -> None:
        return self.plot_seed()

    @deprecated
    def plotInit(self) -> None:
        return self.plot_init()

    @deprecated
    def plotSim(self) -> None:
       return self.plot_sim()
