#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import os
import os.path
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from betse.exceptions import (
    BetseFileException, BetseExceptionSimulation, BetseExceptionParameters)
from betse.science import filehandling as fh
from betse.science.cells import Cells
from betse.science.parameters import Parameters
from betse.science.plot import plot as viz
from betse.science.plot.pipeline import plot_all
from betse.science.sim import Simulator
from betse.science.tissue.handler import TissueHandler
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.science.chemistry.metabolism import MasterOfMetabolism
from betse.science.chemistry.gene import MasterOfGenes


class SimRunner(object):
    '''
    High-level simulation runner encapsulating a single simulation.

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
    def __init__(self, config_filename: str):
        super().__init__()

        # Validate and localize such filename.
        files.die_unless_file(config_filename)
        self._config_filename = config_filename
        self._config_basename = paths.get_basename(self._config_filename)

    def makeWorld(self):
        """
        In order to set up tissue profiles and other geometry-specific features,
        it is necessary to first create and plot the cells data structure. This
        will be loaded into the initialization and simulation runs.

        Parameters
        ----------
        plotWorld : bool, optional
            True if a non-blocking plot of the created cellular world is to be
            displayed immediately after creating such world. Defaults to False,
            in which case no plot will be displayed.
        """

        logs.log_info(
            'Seeding simulation with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename=self._config_filename)     # create an instance of Parameters
        p.I_overlay = False  # force the current overlay to be null
        sim = Simulator(p)   # create an instance of Simulator as it's needed by plotting objects

        cells = Cells(p)  # create an instance of the Cells object
        logs.log_info('Cell cluster is being created...')
        cells.makeWorld(p)  # call function to create the world

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

        if p.sim_eosmosis is True:
            cells.eosmo_tools(p)

        # finish up:

        cells.save_cluster(p)

        logs.log_info('Cell cluster creation complete!')

        if p.turn_all_plots_off is False:
            logs.log_info('Close all plot windows to continue...')
            self.plotWorld()

            plt.show()

        sim.sim_info_report(cells,p)

    def initialize(self):
        '''
        Run an initialization simulation from scratch and save it to the
        initialization cache.
        '''

        logs.log_info(
            'Initializing simulation with configuration file "{}".'.format(
                self._config_basename))

        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.set_time_profile(p.time_profile_init)  # force the time profile to be initialize
        p.run_sim = False # let the simulator know we're just running an initialization

        # cells, _ = fh.loadSim(cells.savedWorld)
        cells = Cells(p)  # create an instance of world

        if files.is_file(cells.savedWorld):
            cells,p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')

            # check to ensure compatibility between original and present sim files:
            self.check_congruency(p_old,p)

        else:
            logs.log_info("Ooops! No such cell cluster file found to load!")

            if p.autoInit is True:
                logs.log_info(
                    'Automatically seeding cell cluster from config file settings...')
                self.makeWorld()  # create an instance of world
                logs.log_info(
                    'Now using cell cluster to run initialization.')
                cells,_ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

            else:
                raise BetseExceptionSimulation(
                    "Run terminated due to missing seed.\n"
                    "Please run 'betse seed' to try again.")

        sim = Simulator(p)   # create an instance of Simulator
        sim.run_sim = False

        # Initialize simulation data structures, run, and save simulation phase
        sim.baseInit_all(cells, p)
        sim.sim_info_report(cells, p)
        sim.run_sim_core(cells, p)

        logs.log_info('Initialization run complete!')
        logs.log_info(
            'The initialization took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        # As colormaps are deleted from the current Parameters instance "p"
        # prior to saving in sim, create a fresh instance of Parameters.
        p = Parameters(config_filename = self._config_filename)

        # Create all enabled plots and animations.
        logs.log_info(
            'When ready, close all of the figure windows to proceed with '
            'scheduled simulation runs.')
        plot_all(cells, sim, p, plot_type='init')

    def simulate(self):
        '''
        Run simulation from a previously saved initialization.
        '''
        logs.log_info(
            'Running simulation with configuration file "{}".'.format(
                self._config_basename))

        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.set_time_profile(p.time_profile_sim)  # force the time profile to be initialize
        p.run_sim = True    # set on the fly a boolean to let simulator know we're running a full simulation
        sim = Simulator(p)   # create an instance of Simulator

        if files.is_file(sim.savedInit):
            sim,cells, p_old = fh.loadSim(sim.savedInit)  # load the initialization from cache
            p.sim_ECM = cells.sim_ECM

            # check to ensure compatibility between original and present sim files:
            self.check_congruency(p_old,p)

        else:
            logs.log_info("No initialization file found to run this simulation!")

            if p.autoInit is True:
                logs.log_info("Automatically running initialization...")
                self.initialize()
                logs.log_info('Now using initialization to run simulation.')
                sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

            elif p.autoInit is False:
                raise BetseExceptionSimulation("Simulation terminated due to missing initialization. Please run "
                                               "an initialization and try again.")

        # Reinitialize save and load directories in case params defines new ones
        # for this sim.
        sim.fileInit(p)

        # Run and save the simulation to the cache.
        sim.sim_info_report(cells, p)

        sim.run_sim_core(cells, p)

        logs.log_info(
            'The simulation took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        # As colormaps are deleted from the current Parameters instance "p"
        # prior to saving in sim, create a fresh instance of Parameters.
        p = Parameters(config_filename=self._config_filename)

        # Create all enabled plots and animations.
        logs.log_info(
            'When ready, close all of the figure windows to end the program.')
        plot_all(cells, sim, p, plot_type='sim')

    def sim_brn(self):
        '''
        Test run a bioenergetics reaction network (without bioelectrics) and
        save it to the initialization cache.
        '''

        logs.log_info(
            'Testing bioenergetics reaction network indicated in configuration file "{}".'.format(
                self._config_basename))

        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters(config_filename=self._config_filename)  # create an instance of Parameters
        p.set_time_profile(p.time_profile_init)  # force the time profile to be initialize
        p.run_sim = False

        # cells, _ = fh.loadSim(cells.savedWorld)
        cells = Cells(p)  # create an instance of world

        if files.is_file(cells.savedWorld):
            cells, p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')

            # check to ensure compatibility between original and present sim files:
            self.check_congruency(p_old, p)

        else:
            logs.log_info("Ooops! No such cell cluster file found to load!")

            if p.autoInit is True:
                logs.log_info(
                    'Automatically seeding cell cluster from config file settings...')
                self.makeWorld()  # create an instance of world
                logs.log_info(
                    'Now using cell cluster to run initialization.')
                cells, _ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

            else:
                raise BetseExceptionSimulation(
                    "Run terminated due to missing seed.\n"
                    "Please run 'betse seed' to try again.")

        sim = Simulator(p)  # create an instance of Simulator

        # Initialize simulation data structures
        sim.baseInit_all(cells, p)

        # create an instance of master of metabolism
        MoM = MasterOfMetabolism(p)

        # initialize it:
        MoM.read_metabo_config(sim, cells, p)

        logs.log_info("Running metabolic reaction network test simulation...")

        MoM.run_core_sim(sim, cells, p)

        logs.log_info(
            'This metabolic network test took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

    def plot_brn(self):

        p = Parameters(config_filename=self._config_filename)  # create an instance of Parameters

        MoM = MasterOfMetabolism(p)

        MoM, cells, _ = fh.loadSim(MoM.savedMoM)

        sim = Simulator(p)  # create an instance of Simulator
        # Initialize simulation data structures
        sim.baseInit_all(cells, p)
        sim.time = MoM.time

        MoM.core.init_saving(cells, p, plot_type = 'init',nested_folder_name='Metabolism')
        MoM.core.export_all_data(sim, cells, p, message = 'for metabolic molecules...')
        MoM.core.plot(sim, cells, p, message = 'for metabolic molecules...')
        MoM.core.anim(sim, cells, p, message = 'for metabolic molecules...')

        if p.turn_all_plots_off is False:
            plt.show()

    def sim_grn(self):

        '''
            Test run a gene regulatory network (without bioelectrics) and save it to the
            initialization cache.
            '''

        logs.log_info(
            'Testing gene regulatory network indicated in configuration file "{}".'.format(
                self._config_basename))

        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters(config_filename=self._config_filename)  # create an instance of Parameters
        p.set_time_profile(p.time_profile_init)  # force the time profile to be initialize
        p.run_sim = False

        # cells, _ = fh.loadSim(cells.savedWorld)
        cells = Cells(p)  # create an instance of world

        if files.is_file(cells.savedWorld):
            cells, p_old = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            logs.log_info('Cell cluster loaded.')

            # check to ensure compatibility between original and present sim files:
            self.check_congruency(p_old, p)

        else:
            logs.log_info("Ooops! No such cell cluster file found to load!")

            if p.autoInit is True:
                logs.log_info(
                    'Automatically seeding cell cluster from config file settings...')
                self.makeWorld()  # create an instance of world
                logs.log_info(
                    'Now using cell cluster to run initialization.')
                cells, _ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

            else:
                raise BetseExceptionSimulation(
                    "Run terminated due to missing seed.\n"
                    "Please run 'betse seed' to try again.")

        sim = Simulator(p)  # create an instance of Simulator

        # Initialize simulation data structures
        sim.baseInit_all(cells, p)

        # create an instance of master of metabolism
        MoG = MasterOfGenes(p)

        # initialize it:
        MoG.read_gene_config(sim, cells, p)

        logs.log_info("Running gene regulatory network test simulation...")

        MoG.run_core_sim(sim, cells, p)

        logs.log_info(
            'This gene regulatory network test took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

    def plot_grn(self):

        p = Parameters(config_filename=self._config_filename)  # create an instance of Parameters

        MoG = MasterOfGenes(p)

        MoG, cells, _ = fh.loadSim(MoG.savedMoG)

        sim = Simulator(p)  # create an instance of Simulator
        # Initialize simulation data structures
        sim.baseInit_all(cells, p)
        sim.time = MoG.time

        MoG.core.init_saving(cells, p, plot_type = 'init',nested_folder_name='GRN')
        MoG.core.export_all_data(sim, cells, p, message = 'for gene products...')
        MoG.core.plot(sim, cells, p, message = 'for gene products...')
        MoG.core.anim(sim, cells, p, message = 'for gene products...')

        if p.turn_all_plots_off is False:
            plt.show()

    def plotInit(self):
        '''
        Load and visualize a previously solved initialization.
        '''
        logs.log_info(
            'Plotting initialization with configuration "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator

        if files.is_file(sim.savedInit):
            sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        else:
            raise BetseExceptionSimulation(
                "Ooops! No such initialization file found to plot!")

        plot_all(cells, sim, p, plot_type='init')

        # run the molecules plots:
        if p.molecules_enabled and sim.molecules is not None:

            # reinit settings for plots, in case they've changed:
            sim.molecules.plot_init(p.molecules_config)

            sim.molecules.init_saving(cells, p, plot_type = 'init')
            sim.molecules.export_all_data(sim, cells, p)
            sim.molecules.plot(sim, cells, p)
            sim.molecules.anim(sim, cells, p)

        if p.metabolism_enabled and sim.metabo is not None:

            sim.metabo.core.init_saving(cells, p, plot_type='init', nested_folder_name='Metabolism')
            sim.metabo.core.export_all_data(sim, cells, p, message = 'for metabolic molecules...')
            sim.metabo.core.plot(sim, cells, p, message = 'for metabolic molecules...')
            sim.metabo.core.anim(sim, cells, p, message = 'for metabolic molecules...')

        if p.grn_enabled and sim.grn is not None:

            sim.grn.core.init_saving(cells, p, plot_type='init', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, p, message = 'for GRN molecules...')
            sim.grn.core.plot(sim, cells, p, message = 'for GRN molecules...')
            sim.grn.core.anim(sim, cells, p, message = 'for GRN molecules...')

        if p.Ca_dyn is True and p.ions_dict['Ca'] == 1:

            sim.endo_retic.init_saving(cells, p, plot_type = 'init', nested_folder_name = 'ER')
            sim.endo_retic.plot_er(sim, cells, p)

        if p.turn_all_plots_off is False:
            plt.show()

    def plotSim(self):
        '''
        Load and visualize a previously solved simulation.
        '''
        logs.log_info(
            'Plotting simulation with configuration "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator

        # If this simulation has yet to be run, fail.
        if not files.is_file(sim.savedSim):
            raise BetseFileException(
                'Simulation cache file "{}" not found to plot '
                '(e.g., due to no simulation having been run).'.format(
                    sim.savedSim))

        # Load the simulation from the cache.
        sim, cells, _ = fh.loadSim(sim.savedSim)
        plot_all(cells, sim, p, plot_type='sim')

        # run the molecules plots:
        if p.molecules_enabled and sim.molecules is not None:

            # reinit settings for plots, in case they've changed:
            sim.molecules.plot_init(p.molecules_config)

            sim.molecules.init_saving(cells, p, plot_type = 'sim')
            sim.molecules.export_all_data(sim, cells, p)
            sim.molecules.plot(sim, cells, p)
            sim.molecules.anim(sim, cells, p)

        if p.metabolism_enabled and sim.metabo is not None:

            sim.metabo.core.init_saving(cells, p, plot_type='sim')
            sim.metabo.core.export_all_data(sim, cells, p)
            sim.metabo.core.plot(sim, cells, p)
            sim.metabo.core.anim(sim, cells, p)

        if p.grn_enabled and sim.grn is not None:

            sim.grn.core.init_saving(cells, p, plot_type='sim', nested_folder_name='GRN')
            sim.grn.core.export_all_data(sim, cells, p, message = 'for GRN molecules...')
            sim.grn.core.plot(sim, cells, p, message = 'for GRN molecules...')
            sim.grn.core.anim(sim, cells, p, message = 'for GRN molecules...')

        if p.Ca_dyn is True and p.ions_dict['Ca'] == 1:

            sim.endo_retic.init_saving(cells, p, plot_type = 'sim', nested_folder_name = 'ER')
            sim.endo_retic.plot_er(sim, cells, p)

        if p.turn_all_plots_off is False:
            plt.show()

    def plotWorld(self):
        '''
        Load and visualize a previously seeded cell cluster.
        '''
        logs.log_info(
            'Plotting cell cluster with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.I_overlay = False # force the current overlay to be false as there's no data for it
        sim = Simulator(p)

        cells = Cells(p)

        if files.is_file(cells.savedWorld):
            cells,_ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            p.sim_ECM = cells.sim_ECM
            logs.log_info('Cell cluster loaded.')
        else:
            raise BetseExceptionSimulation("Ooops! No such cell cluster file found to load!")

        sim.baseInit_all(cells,p)
        dyna = TissueHandler(sim,cells,p)
        dyna.tissueProfiles(sim,cells,p)

        if p.autosave is True:

            images_path = p.init_results
            image_cache_dir = os.path.expanduser(images_path)
            os.makedirs(image_cache_dir, exist_ok=True)
            savedImg = os.path.join(image_cache_dir, 'fig_')

        if p.plot_cell_cluster is True:

            fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(
                p, dyna, cells, clrmap=p.default_cm)

            if p.autosave is True:
                savename10 = savedImg + 'cluster_mosaic' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block = False)

        if p.sim_ECM is True and p.plot_cluster_mask is True:
            plt.figure()
            ax99 = plt.subplot(111)
            plt.imshow(
                np.log10(sim.D_env_weight.reshape(cells.X.shape)),
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()

            cell_edges_flat = p.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(1.0)
            ax99.add_collection(coll)

            plt.title('Logarithm of Environmental Diffusion Weight Matrix')

            if p.autosave is True:
                savename10 = savedImg + 'env_diffusion_weights' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block = False)

            plt.figure()
            plt.imshow(
                cells.maskM,
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()
            plt.title('Cluster Masking Matrix')

            if p.autosave is True:
                savename = savedImg + 'cluster_mask' + '.png'
                plt.savefig(savename,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block = False)


        if p.plot_cell_connectivity is True:
            # plot gj
            fig_x = plt.figure()
            ax_x = plt.subplot(111)

            if p.showCells is True:
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
                plt.show(block = False)

        if p.turn_all_plots_off is False:
            plt.show(block = False)
            plt.show()

        else:
            logs.log_info(
                'Plots exported to init results folder '
                'defined in configuration file "{}".'.format(
                    self._config_basename))

    def check_congruency(self, p_old, p):

        if p_old.config['general options'] != p.config['general options'] or \
           p_old.config['world options'] != p.config['world options']:

            # logs.log_warning('---------------------------------------------------')
            # logs.log_warning('**WARNING!**')
            # logs.log_warning('Important config file options are out of sync ')
            # logs.log_warning('between the seed and this init/sim attempt! ')
            # logs.log_warning('Run "betse seed" again to match the current settings')
            # logs.log_warning(' of this config file.')
            # logs.log_warning('---------------------------------------------------')

                raise BetseExceptionParameters(
                    'Important config file options are out of sync between '
                    'seed and this init/sim attempt! '
                    'Run "betse seed" again to match the current settings of '
                    'this config file.')

