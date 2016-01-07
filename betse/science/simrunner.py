#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


import os, os.path, time
import matplotlib.pyplot as plt
import numpy as np
from betse.exceptions import BetseExceptionSimulation, BetseExceptionParameters
from betse.science import filehandling as fh
from betse.science import visualize as viz
from betse.science.animation.animation import (
    AnimateCellData,
    AnimateGJData,
    AnimateCurrent,
    AnimateField,
    AnimateVelocity,
    AnimateDeformation,
    AnimateEnv,
    AnimateMem,
    AnimateDyeData,
)
from betse.science.cells import Cells
from betse.science.parameters import Parameters
from betse.science.sim import Simulator
from betse.science.tissue.handler import TissueHandler
from betse.util.io import loggers
from betse.util.path import files, paths
from matplotlib.collections import LineCollection, PolyCollection


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

        loggers.log_info(
            'Seeding simulation with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename=self._config_filename)     # create an instance of Parameters
        p.I_overlay = False  # force the current overlay to be null
        sim = Simulator(p)   # create an instance of Simulator as it's needed by plotting objects

        if p.sim_ECM is False:

            cells = Cells(p,worldtype='basic')  # create an instance of world
            cells.containsECM = False
            loggers.log_info('Cell cluster is being created...')
            cells.makeWorld(p)     # call function to create the world

            # define the tissue and boundary profiles for plotting:
            loggers.log_info('Defining tissue and boundary profiles...')
            sim.baseInit(cells,p)
            dyna = TissueHandler(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)

            cells.redo_gj(dyna,p)  # redo gap junctions to isolate different tissue types

            if p.fluid_flow is True or p.deformation is True: # if user desires fluid flow:

                # make a laplacian and solver for discrete transfers on closed, irregular cell network
                loggers.log_info('Creating cell network Poisson solver...')
                cells.graphLaplacian(p)

            cells.save_cluster(p)

            loggers.log_info('Cell cluster creation complete!')

            if p.turn_all_plots_off is False:
                loggers.log_info('Close all plot windows to continue...')
                self.plotWorld()

        else:

            cells = Cells(p,worldtype='full')  # create an instance of world
            cells.containsECM = True
            loggers.log_info('Cell cluster is being created...')
            cells.makeWorld(p)     # call function to create the world

            # define the tissue and boundary profiles for plotting:
            loggers.log_info('Defining tissue and boundary profiles...')
            sim.baseInit_ECM(cells,p)
            dyna = TissueHandler(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)

            cells.redo_gj(dyna,p)  # redo gap junctions to isolate different tissue types

            # make a laplacian and solver for discrete transfers on closed, irregular cell network
            if p.fluid_flow is True or p.deformation is True:
                # loggers.log_info('Creating cell network Poisson solvers...')
                cells.graphLaplacian(p)

            cells.save_cluster(p)

            loggers.log_info('Cell cluster creation complete!')

            if p.turn_all_plots_off is False:
                loggers.log_info('Close all plot windows to continue...')
                self.plotWorld()

                plt.show()

    def initialize(self):
        '''
        Run an initialization simulation from scratch and save it to the
        initialization cache.
        '''

        loggers.log_info(
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
            loggers.log_info('Cell cluster loaded.')

            #FIXME: This same test is duplicate in simulate(). Diadem trinkets!
            if p_old.config['general options'] != p.config['general options'] or \
               p_old.config['world options'] != p.config['world options'] or \
               p_old.config['tissue profile definition'] != p.config['tissue profile definition']:
                raise BetseExceptionParameters(
                    'Important config file options are out of sync between seed and this init attempt! '
                    'Run "betse seed" again to match the current settings of this config file.')

        else:
            loggers.log_info("Ooops! No such cell cluster file found to load!")

            if p.autoInit is True:
                loggers.log_info("Automatically seeding cell cluster from config file settings...")
                self.makeWorld()  # create an instance of world
                loggers.log_info('Now using cell cluster to run initialization.')
                cells,_ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache

            else:
                raise BetseExceptionSimulation(
                    "Run terminated due to missing seed.\n"
                    "Please run 'betse seed' to try again.")

        sim = Simulator(p)   # create an instance of Simulator
        sim.run_sim = False

        if p.sim_ECM is False:
            sim.baseInit(cells, p)   # initialize simulation data structures
            # sim.tissueInit(cells,p)
            sim.runSim(cells,p)     # run and save the initialization

        elif p.sim_ECM is True:
            sim.baseInit_ECM(cells, p)   # initialize simulation data structures
            # sim.tissueInit(cells,p)
            sim.runSim_ECM(cells,p)     # run and save the initialization

        loggers.log_info('Initialization run complete!')
        loggers.log_info(
            'The initialization took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        if p.turn_all_plots_off is False:
            # as colormaps are deleted from p prior to saving in sim, create a fresh instance of Parameters:
            p = Parameters(config_filename = self._config_filename)

            loggers.log_info(
                'When ready, close all of the figure windows to proceed with '
                'scheduled simulation runs.')

            plot_sim(cells, sim, p, plot_type='init')
            plt.show()


    def simulate(self):
        '''
        Run simulation from a previously saved initialization.
        '''
        loggers.log_info(
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

            #FIXME: This same test is duplicate in initialize(). Emerald sky!
            if p_old.config['general options'] != p.config['general options'] or \
               p_old.config['world options'] != p.config['world options'] or \
               p_old.config['tissue profile definition'] != p.config['tissue profile definition']:
                raise BetseExceptionParameters(
                    'Important config file options are out of sync between the seed and this sim attempt! '
                    'Run "betse seed" and "betse init" again to match the current settings of this config file.')

        else:
            loggers.log_info("No initialization file found to run this simulation!")

            if p.autoInit is True:
                loggers.log_info("Automatically running initialization...")
                self.initialize()
                loggers.log_info('Now using initialization to run simulation.')
                sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

            elif p.autoInit is False:
                raise BetseExceptionSimulation("Simulation terminated due to missing initialization. Please run"
                                               "an initialization and try again.")

        sim.fileInit(p)   # reinitialize save and load directories in case params defines new ones for this sim

        # Run and optionally save the simulation to the cache.
        if p.sim_ECM is False:
            sim.runSim(cells,p,save=True)
        else:
            sim.runSim_ECM(cells,p,save=True)

        loggers.log_info(
            'The simulation took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))
        loggers.log_info(
            'When ready, close all of the figure windows to end the program.')

        if p.turn_all_plots_off is False:
            #FIXME: Is there no better way? As "Parameters" is a fairly heavy
            #class, it'd be great to reuse the existing "p" object. Perhaps we
            #could keep around colormaps somehow? Infants raised in the jungle!

            # As colormaps are deleted from p prior to saving in sim, create a
            # fresh instance of Parameters.
            p = Parameters(config_filename=self._config_filename)
            plot_sim(cells, sim, p, plot_type='sim')
            plt.show()


    def plotInit(self):
        '''
        Load and visualize a previously solved initialization.
        '''
        loggers.log_info(
            'Plotting initialization with configuration "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator

        if files.is_file(sim.savedInit):
            sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        else:
            raise BetseExceptionSimulation(
                "Ooops! No such initialization file found to plot!")

        plot_sim(cells, sim, p, plot_type='init')

        if p.turn_all_plots_off is False:
            plt.show()

    def plotSim(self):
        '''
        Load and visualize a previously solved simulation.
        '''
        loggers.log_info(
            'Plotting simulation with configuration "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator

        if files.is_file(sim.savedSim):
            sim,cells,_ = fh.loadSim(sim.savedSim)  # load the simulation from cache
        else:
            raise BetseExceptionSimulation(
                "Ooops! No such simulation file found to plot!")

        plot_sim(cells, sim, p, plot_type='sim')

        if p.turn_all_plots_off is False:
            plt.show()

    def plotWorld(self):

        loggers.log_info(
            'Plotting cell cluster with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.I_overlay = False # force the current overlay to be false as there's no data for it
        sim = Simulator(p)

        cells = Cells(p,worldtype='basic')

        if files.is_file(cells.savedWorld):
            cells,_ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            p.sim_ECM = cells.sim_ECM
            loggers.log_info('Cell cluster loaded.')
        else:
            raise BetseExceptionSimulation("Ooops! No such cell cluster file found to load!")

        if p.sim_ECM is False:
            sim.baseInit(cells,p)
            dyna = TissueHandler(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)
        else:
            sim.baseInit_ECM(cells,p)
            dyna = TissueHandler(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)

        if p.autosave is True:

            images_path = p.init_results
            image_cache_dir = os.path.expanduser(images_path)
            os.makedirs(image_cache_dir, exist_ok=True)
            savedImg = os.path.join(image_cache_dir, 'fig_')

        fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(
            p, dyna, cells, clrmap=p.default_cm)

        if p.autosave is True:
            savename10 = savedImg + 'cluster_mosaic' + '.png'
            plt.savefig(savename10,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block = False)

        if p.sim_ECM is True:
            plt.figure()
            ax99 = plt.subplot(111)
            plt.imshow(
                np.log10(sim.D_env_weight.reshape(cells.X.shape)),
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()

            cell_edges_flat = cells.um*cells.mem_edges_flat
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

            if p.turn_all_plots_off is False:
                plt.show(block = False)

        # plot gj
        fig_x = plt.figure()
        ax_x = plt.subplot(111)

        if p.showCells is True:
            base_points = np.multiply(cells.cell_verts, p.um)
            col_cells = PolyCollection(base_points, facecolors='k', edgecolors='none')
            col_cells.set_alpha(0.3)
            ax_x.add_collection(col_cells)

            # cell_edges_flat = p.um*cells.mem_edges_flat
            # coll = LineCollection(cell_edges_flat,color='k',linewidth=0.5)
            # coll.set_alpha(0.5)
            # ax_x.add_collection(coll)

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
            plt.show()

        else:
            loggers.log_info(
                'Plots exported to init results folder '
                'defined in configuration file "{}".'.format(
                    self._config_basename))

#FIXME: This single function appears to constitute over two thirds of this
#entire module! It's also one of the more intense functions in any language
#I've ever seen. That's not necessarily a bad thing; it's pretty cool in its
#complexity, and I can't believe it actually works. That said, it doesn't seem
#to have much to do with high-level simulation running. It's more of a
#low-level plots manager, right? Perhaps we could shove this away into a
#"betse.science.plot" module. Shivering reindeer pawing at the door!
def plot_sim(cells, sim, p, plot_type: str = 'init'):

    if p.autosave is True:
        if plot_type == 'sim':
            images_path = p.sim_results
            p.plot_type = 'sim'

        elif plot_type == 'init':
            images_path = p.init_results
            p.plot_type = 'init'

        image_cache_dir = os.path.expanduser(images_path)
        os.makedirs(image_cache_dir, exist_ok=True)
        savedImg = os.path.join(image_cache_dir, 'fig_')

    # check that the plot cell is in range of the available cell indices:
    if p.plot_cell not in cells.cell_i:
        raise BetseExceptionParameters(
            'The "plot cell" defined in the "results" section of your '
            'configuration file does not exist in your cluster. '
            'Choose a plot cell number smaller than the maximum cell number.')

    if p.sim_ECM is True:
        plot_cell_ecm = cells.cell_to_mems[p.plot_cell][0]  # convert from cell to mem index

    else:
        plot_cell_ecm = p.plot_cell

    if p.exportData is True:
        viz.exportData(cells, sim, p)

    #-------------------------------------------------------------------------------------------------------------------
    #               SINGLE CELL DATA GRAPHS
    #-------------------------------------------------------------------------------------------------------------------

    # plot-cell sodium concentration vs time:

    if p.plot_single_cell_graphs is True:
        figConcsNa, axConcsNa = viz.plotSingleCellCData(
            sim.cc_time, sim.time, sim.iNa, p.plot_cell, fig=None,
            ax=None, lncolor='g', ionname='Na+')

        titNa = 'Sodium concentration in cell ' + str(p.plot_cell)
        axConcsNa.set_title(titNa)

        if p.autosave is True:
            savename1 = savedImg + 'concNa_time' + '.png'
            plt.savefig(savename1, dpi=300, format='png', transparent=True)

    # plot-cell potassium concentration vs time:

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        figConcsK, axConcsK = viz.plotSingleCellCData(
            sim.cc_time, sim.time, sim.iK, p.plot_cell, fig=None,
            ax=None, lncolor='b', ionname='K+')

        titK = 'Potassium concentration in cell ' + str(p.plot_cell)
        axConcsK.set_title(titK)

        if p.autosave is True:
            savename1 = savedImg + 'concK_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot-cell anion (bicarbonate) concentration vs time:

        figConcsM, axConcsM = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,p.plot_cell,fig=None,
             ax=None,lncolor='r',ionname='M-')

        titM = 'M Anion concentration in cell ' + str(p.plot_cell)
        axConcsM.set_title(titM)

        if p.autosave is True:
            savename1 = savedImg + 'concM_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot single cell Vmem vs time:

        figVt, axVt = viz.plotSingleCellVData(sim,plot_cell_ecm,p,fig=None,ax=None,lncolor='k')
        titV = 'Voltage (Vmem) in cell ' + str(p.plot_cell)
        axVt.set_title(titV)

        if p.autosave is True:
            savename2 = savedImg + 'Vmem_time' + '.png'
            plt.savefig(savename2,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot fast-Fourier-transform (fft) of Vmem:

        figFFT, axFFT = viz.plotFFT(sim.time,sim.vm_time,plot_cell_ecm,lab="Power")
        titFFT = 'Fourier transform of Vmem in cell ' + str(p.plot_cell)
        axFFT.set_title(titFFT)

        if p.autosave is True:
            savename = savedImg + 'FFT_time' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot rate of Na-K-ATPase pump vs time:
        figNaK = plt.figure()
        axNaK = plt.subplot()

        if p.sim_ECM is False:
            pump_rate = [pump_array[p.plot_cell] for pump_array in sim.rate_NaKATP_time]

        else:
            pump_rate = [pump_array[plot_cell_ecm] for pump_array in sim.rate_NaKATP_time]

        axNaK.plot(sim.time,pump_rate)

        axNaK.set_xlabel('Time [s]')
        axNaK.set_ylabel('Pumping Rate of Na+ Out of Cell [mol/(m2 s)]')

        axNaK.set_title('Rate of NaK-ATPase pump in cell: ' + str(p.plot_cell) )

        if p.autosave is True:
            savename = savedImg + 'NaKATPaseRaTE_' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        #--------------------------------------------------------

        # plot-cell trans-membrane current vs time:
        figI = plt.figure()
        axI = plt.subplot(111)

        if p.sim_ECM is False:
            Imem = [memArray[p.plot_cell] for memArray in sim.I_mem_time]

        else:
            Imem = []  # initialize a total cell current storage vector
            mems_for_plotcell = cells.cell_to_mems[p.plot_cell]  # get membranes for the plot cell

            for t in range(len(sim.time)):
                memArray = sim.I_mem_time[t]
                # get the current components at each membrane (net current, not density)
                Ixo = memArray[mems_for_plotcell]*cells.mem_vects_flat[mems_for_plotcell,2]*cells.mem_sa[mems_for_plotcell]
                Iyo = memArray[mems_for_plotcell]*cells.mem_vects_flat[mems_for_plotcell,3]*cells.mem_sa[mems_for_plotcell]
                # add components of current at each membrane (this takes account for in-out directionality)
                Ix = np.sum(Ixo)
                Iy = np.sum(Iyo)

                # get the total magnitude of net current and divide by cell surface area to return to density:
                Io = np.sqrt(Ix**2 + Iy**2)/cells.cell_sa[p.plot_cell]
                Imem.append(Io)

        axI.plot(sim.time,Imem)
        axI.set_title('Transmembrane current density for cell ' + str(p.plot_cell) )
        axI.set_xlabel('Time [s]')
        axI.set_ylabel('Current density [A/m2]')

        if p.autosave is True:
            savename = savedImg + 'Imem_time' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # optional 1D plots--------------------------------------------------------------------------------------------

        # hydrostatic pressure in cells:

        if p.deform_osmo is True:
            p_hydro = [arr[p.plot_cell] for arr in sim.P_cells_time]
            figOP = plt.figure()
            axOP = plt.subplot(111)
            axOP.plot(sim.time,p_hydro)
            axOP.set_xlabel('Time [s]')
            axOP.set_ylabel('Hydrostatic Pressure [Pa]')
            axOP.set_title('Hydrostatic pressure in cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'HydrostaticP_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        # plot-cell calcium vs time (if Ca enabled in ion profiles):

        if p.ions_dict['Ca'] ==1:
            figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,p.plot_cell,fig=None,
                 ax=None,lncolor='g',ionname='Ca2+ cell')
            titCa =  'Cytosolic Ca2+ in cell index ' + str(p.plot_cell)
            axA.set_title(titCa)

            if p.autosave is True:
                savename3 = savedImg + 'cytosol_Ca_time' + '.png'
                plt.savefig(savename3,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            if p.Ca_dyn == 1:
                figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,p.plot_cell,fig=None,
                     ax=None,lncolor='b',ionname='Ca2+ cell')
                titER =  'ER Ca2+ in cell index ' + str(p.plot_cell)
                axD.set_title(titER)

                if p.autosave is True:
                    savename4 = savedImg + 'ER_Ca_time' + '.png'
                    plt.savefig(savename4,dpi=300,format='png',transparent=True)

                if p.turn_all_plots_off is False:
                    plt.show(block=False)

                figPro, axPro = viz.plotSingleCellData(sim.time, sim.cIP3_time,p.plot_cell, lab='IP3 [mmol/L]')
                titIP3 =  'IP3 in cell index ' + str(p.plot_cell)
                axPro.set_title(titIP3)

                if p.autosave is True:
                    savename5 = savedImg + 'IP3_time' + '.png'
                    plt.savefig(savename5,dpi=300,format='png',transparent=True)

                if p.turn_all_plots_off is False:
                    plt.show(block=False)

        # osmotic and/or electrostatic pressure in cell
        if p.deform_electro is True:
            f_electro = [arr[p.plot_cell] for arr in sim.P_electro_time]
            figPE = plt.figure()
            axPE = plt.subplot(111)
            axPE.plot(sim.time,f_electro)
            axPE.set_xlabel('Time [s]')
            axPE.set_ylabel('Electrostatic Pressure [Pa]')
            axPE.set_title('Electrostatic pressure in cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'ElectrostaticP_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        if p.deform_osmo is True:
            p_osmo = [arr[p.plot_cell] for arr in sim.osmo_P_delta_time]
            figOP = plt.figure()
            axOP = plt.subplot(111)
            axOP.plot(sim.time,p_osmo)
            axOP.set_xlabel('Time [s]')
            axOP.set_ylabel('Osmotic Pressure [Pa]')
            axOP.set_title('Osmotic pressure in cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'OsmoticP_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        # total displacement in cell
        if p.deformation is True and sim.run_sim is True:
            # extract time-series deformation data for the plot cell:
            dx = np.asarray([arr[p.plot_cell] for arr in sim.dx_cell_time])
            dy = np.asarray([arr[p.plot_cell] for arr in sim.dy_cell_time])

            # get the total magnitude:
            disp = np.sqrt(dx**2 + dy**2)

            figD = plt.figure()
            axD = plt.subplot(111)
            axD.plot(sim.time,p.um*disp)
            axD.set_xlabel('Time [s]')
            axD.set_ylabel('Displacement [um]')
            axD.set_title('Displacement of cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'Displacement_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------
    #                       2D Data Map Plotting
    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_venv is True and p.sim_ECM is True:
        plt.figure()
        venv_plt = plt.imshow(
            1000*sim.v_env.reshape(cells.X.shape),
            origin='lower',
            extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
            cmap=p.default_cm)

        if p.autoscale_venv is False:
            venv_plt.set_clim(p.venv_min_clr, p.venv_max_clr)

        plt.colorbar()
        plt.title('Environmental Voltage [mV]')

        if p.autosave is True:
            savename10 = savedImg + 'Final_environmental_V' + '.png'
            plt.savefig(savename10,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.plot_rho2d is True:
        if p.sim_ECM is True:
            plt.figure()
            plt.imshow(
                sim.rho_env.reshape(cells.X.shape)*(1/p.ff_env),
                origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()
            plt.title('Environmental Charge Density [C/m3]')

            if p.autosave is True:
                savename10 = savedImg + 'Final_environmental_charge' + '.png'
                plt.savefig(savename10,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        figX, axX, cbX = viz.plotPolyData(
            sim, cells, p,
            zdata=(sim.rho_cells)*(1/p.ff_cell),
            number_cells=p.enumerate_cells,
            clrAutoscale=p.autoscale_rho,
            clrMin=p.rho_min_clr,
            clrMax=p.rho_max_clr,
            clrmap=p.default_cm,
            current_overlay=p.I_overlay,
            plotIecm=p.IecmPlot,
        )

        figX.suptitle('Final Cell Charge Density',fontsize=14, fontweight='bold')
        axX.set_xlabel('Spatial distance [um]')
        axX.set_ylabel('Spatial distance [um]')
        cbX.set_label('Net Charge Density [C/m3]')

        if p.autosave is True:
            savename9 = savedImg + 'final_cellCharge' + '.png'
            plt.savefig(savename9,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.plot_vcell2d is True and p.sim_ECM is True:
        figX, axX, cbX = viz.plotPolyData(
            sim, cells, p,
            zdata=sim.vcell_time[-1]*1e3,
            number_cells=p.enumerate_cells,
            clrAutoscale=p.autoscale_vcell,
            clrMin=p.vcell_min_clr,
            clrMax=p.vcell_max_clr,
            clrmap=p.default_cm,
            current_overlay=p.I_overlay,
            plotIecm=p.IecmPlot,
        )

        figX.suptitle('Final Cell Voltage',fontsize=14, fontweight='bold')
        axX.set_xlabel('Spatial distance [um]')
        axX.set_ylabel('Spatial distance [um]')
        cbX.set_label('Voltage mV')

        if p.autosave is True:
            savename9 = savedImg + 'final_cellVoltage' + '.png'
            plt.savefig(savename9,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_vm2d is True:
        if p.sim_ECM is True:
            figV, axV, cbV = viz.plotHetMem(
                sim, cells, p,
                zdata=1000*sim.vm_Matrix[-1],
                number_cells=p.enumerate_cells,
                clrAutoscale=p.autoscale_Vmem,
                clrMin=p.Vmem_min_clr,
                clrMax=p.Vmem_max_clr,
                clrmap=p.default_cm,
                edgeOverlay=p.showCells,
                number_ecm=p.enumerate_cells,
                current_overlay=p.I_overlay,
                plotIecm=p.IecmPlot,
            )

        else:
            figV, axV, cbV = viz.plotPolyData(
                sim, cells, p,
                zdata=1000*sim.vm_time[-1],
                clrAutoscale=p.autoscale_Vmem,
                clrMin=p.Vmem_min_clr,
                clrMax=p.Vmem_max_clr,
                number_cells=p.enumerate_cells,
                clrmap=p.default_cm,
                current_overlay=p.I_overlay,
                plotIecm=p.IecmPlot,
            )

        figV.suptitle('Final Vmem',fontsize=14, fontweight='bold')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        if p.autosave is True:
            savename5 = savedImg + 'final_Vmem_2D' + '.png'
            plt.savefig(savename5,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.GHK_calc is True:
        figV_ghk, axV_ghk, cbV_ghk = viz.plotPolyData(
            sim, cells, p,
            zdata=1000*sim.vm_GHK_time[-1],
            clrAutoscale=p.autoscale_Vmem,
            clrMin=p.Vmem_min_clr,
            clrMax=p.Vmem_max_clr,
            number_cells=p.enumerate_cells,
            clrmap=p.default_cm,
            current_overlay=p.I_overlay,
            plotIecm=p.IecmPlot,
        )

        figV_ghk.suptitle('Final Vmem using Goldman Equation',fontsize=14, fontweight='bold')
        axV_ghk.set_xlabel('Spatial distance [um]')
        axV_ghk.set_ylabel('Spatial distance [um]')
        cbV_ghk.set_label('Voltage [mV]')

        if p.autosave is True:
            savename5 = savedImg + 'final_Vmem_GHK_2D' + '.png'
            plt.savefig(savename5,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------
    if p.plot_ip32d is True and p.scheduled_options['IP3'] != 0 or \
       p.Ca_dyn != 0:
        if p.sim_ECM is False:
            figIP3, axIP3, cbIP3 = viz.plotPolyData(
                sim, cells, p,
                zdata=sim.cIP3_time[-1]*1e3,
                number_cells=p.enumerate_cells,
                clrAutoscale=p.autoscale_IP3,
                clrMin=p.IP3_min_clr,
                clrMax=p.IP3_max_clr,
                clrmap=p.default_cm,
            )

        else:
            figIP3 = plt.figure()
            axIP3 = plt.subplot(111)

            ip3Env = sim.cIP3_env*1e3
            ip3Cell = sim.cIP3*1e3

            bkgPlot = axIP3.imshow(ip3Env.reshape(cells.X.shape),origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,
                p.um*cells.ymax],cmap=p.default_cm)

            points = np.multiply(cells.cell_verts, p.um)

            coll = PolyCollection(points, array=ip3Cell, cmap=p.default_cm, edgecolors='none')
            axIP3.add_collection(coll)
            axIP3.axis('equal')

            # Add a colorbar for the PolyCollection
            maxvala = np.max(ip3Cell,axis=0)
            maxvalb = np.max(ip3Env,axis=0)
            minvala = np.min(ip3Cell,axis=0)
            minvalb = np.min(ip3Env,axis=0)

            #FIXME: Consider using Python's builtin min() and max() functions:
            #    maxval = max(maxvala, maxvalb)
            #    minval = min(minvala, minvalb)
            #Sunrise over a descending moon!

            if maxvala > maxvalb:
                maxval = maxvala
            else:
                maxval = maxvalb

            if minvala < minvalb:
                minval = minvala
            else:
                minval = minvalb

            if p.autoscale_IP3 is True:
                coll.set_clim(minval,maxval)
                bkgPlot.set_clim(minval,maxval)
                cbIP3 = figIP3.colorbar(coll)

            else:
                coll.set_clim(p.IP3_min_clr,p.IP3_max_clr)
                bkgPlot.set_clim(p.IP3_min_clr,p.IP3_max_clr)
                cbIP3 = figIP3.colorbar(coll)

            xmin = cells.xmin*p.um
            xmax = cells.xmax*p.um
            ymin = cells.ymin*p.um
            ymax = cells.ymax*p.um

            axIP3.axis([xmin,xmax,ymin,ymax])

        axIP3.set_title('Final IP3 concentration')
        axIP3.set_xlabel('Spatial distance [um]')
        axIP3.set_ylabel('Spatial distance [um]')
        cbIP3.set_label('Concentration umol/L')

        if p.autosave is True:
            savename6 = savedImg + 'final_IP3_2D' + '.png'
            plt.savefig(savename6,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_dye2d is True and p.voltage_dye == 1:
        if p.sim_ECM is False:

            figVdye, axVdye, cbVdye = viz.plotPolyData(sim, cells,p,zdata=sim.cDye_time[-1]*1e3,
                number_cells=p.enumerate_cells,clrAutoscale = p.autoscale_Dye,
                clrMin = p.Dye_min_clr, clrMax = p.Dye_max_clr, clrmap = p.default_cm)


        else:
            figVdye = plt.figure()
            axVdye = plt.subplot(111)

            dyeEnv = sim.cDye_env*1e3
            dyeCell = sim.cDye_cell*1e3

            bkgPlot = axVdye.imshow(dyeEnv.reshape(cells.X.shape),origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,
                p.um*cells.ymax],cmap=p.default_cm)

            points = np.multiply(cells.cell_verts, p.um)

            coll = PolyCollection(points, array=dyeCell, cmap=p.default_cm, edgecolors='none')
            axVdye.add_collection(coll)
            axVdye.axis('equal')

            # Add a colorbar for the PolyCollection
            maxvala = np.max(dyeCell,axis=0)
            maxvalb = np.max(dyeEnv,axis=0)
            minvala = np.min(dyeCell,axis=0)
            minvalb = np.min(dyeEnv,axis=0)

            if maxvala > maxvalb:
                maxval = maxvala
            else:
                maxval = maxvalb

            if minvala < minvalb:
                minval = minvala
            else:
                minval = minvalb

            if p.autoscale_Dye is True:
                coll.set_clim(minval,maxval)
                bkgPlot.set_clim(minval,maxval)
                cbVdye = figVdye.colorbar(coll)

            elif p.autoscale_Dye is False:
                coll.set_clim(p.Dye_min_clr,p.Dye_max_clr)
                bkgPlot.set_clim(p.Dye_min_clr,p.Dye_max_clr)
                cbVdye = figVdye.colorbar(coll)

            xmin = cells.xmin*p.um
            xmax = cells.xmax*p.um
            ymin = cells.ymin*p.um
            ymax = cells.ymax*p.um

            axVdye.axis([xmin,xmax,ymin,ymax])

        axVdye.set_title('Final Morphogen Concentration')
        axVdye.set_xlabel('Spatial distance [um]')
        axVdye.set_ylabel('Spatial distance [um]')
        cbVdye.set_label('Concentration umol/L')

        if p.autosave is True:
            savename7 = savedImg + 'final_morphogen_2D' + '.png'
            plt.savefig(savename7,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_ca2d is True and p.ions_dict['Ca'] == 1:
        figCa, axCa, cbCa = viz.plotPolyData(sim,cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,
            number_cells= p.enumerate_cells, clrAutoscale = p.autoscale_Ca,
            clrMin = p.Ca_min_clr, clrMax = p.Ca_max_clr, clrmap = p.default_cm)

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if p.autosave is True:
            savename8 = savedImg + 'final_Ca_2D' + '.png'
            plt.savefig(savename8,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.plot_pH2d is True and p.ions_dict['H'] == 1:
        pHdata = -np.log10(sim.cc_time[-1][sim.iH])

        figH, axH, cbH = viz.plotPolyData(sim,cells,p,zdata=pHdata,
            number_cells= p.enumerate_cells, clrAutoscale = p.autoscale_pH,
            clrMin = p.pH_min_clr, clrMax = p.pH_max_clr, clrmap = p.default_cm)

        axH.set_title('Final cytosolic pH')
        axH.set_xlabel('Spatial distance [um]')
        axH.set_ylabel('Spatial distance [um]')
        cbH.set_label('pH')

        if p.autosave is True:
            savename8 = savedImg + 'final_pH_2D' + '.png'
            plt.savefig(savename8,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_I2d is True:
        figI, axI, cbI = viz.streamingCurrent(
            sim, cells, p,
            plot_Iecm=False,
            clrAutoscale=p.autoscale_I2d,
            clrMin=p.I_min_clr,
            clrMax=p.I_max_clr,
            clrmap=p.background_cm,
            edgeOverlay=p.showCells,
            number_cells=p.enumerate_cells,
        )

        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current Density [A/m2]')

        if p.autosave is True:
            savename10 = savedImg + 'Final_Current_gj' + '.png'
            plt.savefig(savename10,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if p.sim_ECM is True:
            figI2, axI2, cbI2 = viz.streamingCurrent(
                sim, cells, p,
                plot_Iecm=True,
                clrAutoscale=p.autoscale_I2d,
                clrMin=p.I_min_clr,
                clrMax=p.I_max_clr,
                clrmap=p.background_cm,
                edgeOverlay=p.showCells,
                number_cells=p.enumerate_cells,
            )

            axI2.set_xlabel('Spatial distance [um]')
            axI2.set_ylabel('Spatial distance [um]')
            cbI2.set_label('Current Density [A/m2]')

            if p.autosave is True:
                savename11 = savedImg + 'Final_Current_extracellular' + '.png'
                plt.savefig(savename11,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_Efield is True:

        if p.sim_ECM is True:

            viz.plotVectField(sim.E_env_x,sim.E_env_y,cells,p,plot_ecm = True,
                title='Final Electric Field', cb_title = 'Electric Field [V/m]',
                colorAutoscale = p.autoscale_Efield, minColor = p.Efield_min_clr,
                maxColor = p.Efield_max_clr)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            if p.autosave is True:
                if p.sim_ECM is True:
                    savename = savedImg + 'Final_Electric_Field_ECM' + '.png'
                    plt.savefig(savename,format='png',transparent=True)

        viz.plotVectField(
            sim.E_gj_x,sim.E_gj_y,cells,p,plot_ecm = False,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=p.autoscale_Efield,
            minColor=p.Efield_min_clr,
            maxColor=p.Efield_max_clr,
        )

        if p.autosave is True:
            savename = savedImg + 'Final_Electric_Field_GJ' + '.png'
            plt.savefig(savename,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_P is True and p.deform_osmo is True:

        figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=sim.P_cells,number_cells=p.enumerate_cells,
        clrAutoscale = p.autoscale_P, clrMin = p.P_min_clr, clrMax = p.P_max_clr, clrmap = p.default_cm)

        axP.set_title('Final Hydrostatic Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure [Pa]')

        if p.autosave is True:
            savename13 = savedImg + 'final_P_2D_gj' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_osmoP is True and p.deform_osmo is True:

        osmo_P = sim.osmo_P_delta

        figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=osmo_P,number_cells=p.enumerate_cells,
        clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr,
            clrmap = p.default_cm)

        axP.set_title('Final Osmotic Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure Difference Cell Interior vs Exterior [Pa]')

        if p.autosave is True:
            savename13 = savedImg + 'final_osmoP_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.plot_osmoP is True and p.deform_electro is True:
        osmo_P = sim.P_electro

        figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=osmo_P,number_cells=p.enumerate_cells,
        clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr,
            clrmap = p.default_cm)

        axP.set_title('Final Electrostatic Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure Difference Cell Interior vs Exterior [Pa]')

        if p.autosave is True:
            savename13 = savedImg + 'final_electroP_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #---- Forces ------------------------------------------------------------

    if p.deform_electro is True:
        viz.plotVectField(
            (1/p.um)*sim.F_electro_x,
            (1/p.um)*sim.F_electro_y,
            cells, p,
            plot_ecm=False,
            title='Final Electrostatic Body Force',
            cb_title='Body Force [N/cm3]',
            colorAutoscale=p.autoscale_force,
            minColor=p.force_min_clr,
            maxColor=p.force_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_electroF_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.deform_osmo is True:
        viz.plotVectField(
            (1/p.um)*sim.F_hydro_x,
            (1/p.um)*sim.F_hydro_y,
            cells, p,
            plot_ecm = False,
            title='Final Hydrostatic Pressure Induced Body Force',
            cb_title='Body Force [N/cm3]',
            colorAutoscale=p.autoscale_force,
            minColor=p.force_min_clr,
            maxColor=p.force_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_hydroF_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.deformation is True and sim.run_sim is True:
        viz.plotStreamField(
            p.um*sim.dx_cell_time[-1],
            p.um*sim.dy_cell_time[-1],
            cells, p,
            plot_ecm=False,
            title='Final Displacement of Cell Collective',
            cb_title='Displacement [um]',
            show_cells=p.showCells,
            colorAutoscale=p.autoscale_Deformation_ani,
            minColor=p.Deformation_ani_min_clr,
            maxColor=p.Deformation_ani_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_displacement_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.plot_Vel is True and p.fluid_flow is True and p.deform_electro is True:
        viz.plotStreamField(
            (1e9)*sim.u_cells_x,
            (1e9)*sim.u_cells_y,
            cells, p,
            plot_ecm=False,
            title='Final Fluid Velocity in Cell Collective',
            cb_title='Velocity [nm/s]',
            colorAutoscale=p.autoscale_Vel,
            minColor=p.Vel_min_clr,
            maxColor=p.Vel_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_vel_2D_gj' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if p.sim_ECM is True:
            viz.plotStreamField((1e9)*sim.u_env_x,(1e9)*sim.u_env_y,cells,p,plot_ecm = True,
                title='Final Fluid Velocity in Cell Collective', cb_title = 'Velocity [nm/s]',
                colorAutoscale = p.autoscale_Vel,minColor = p.Vel_min_clr,maxColor = p.Vel_max_clr)

            if p.autosave is True:
                savename13 = savedImg + 'final_vel_2D_env' + '.png'
                plt.savefig(savename13,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    if p.gj_flux_sensitive is True or p.v_sensitive_gj is True:
        # viz.plotMemData(cells,p,zdata=sim.rho_gj,clrmap=p.default_cm)
        fig_x = plt.figure()
        ax_x = plt.subplot(111)
        con_segs = cells.nn_edges
        connects = p.um*np.asarray(con_segs)
        collection = LineCollection(connects, array=sim.gjopen, cmap= p.default_cm, linewidths=2.0)
        ax_x.add_collection(collection)
        cb = fig_x.colorbar(collection)
        plt.axis('equal')
        plt.axis([cells.xmin*p.um,cells.xmax*p.um,cells.ymin*p.um,cells.ymax*p.um])

        cb.set_label('Relative Permeability')
        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Final Gap Junction Relative Permeability')

        if p.autosave is True:
            savename = savedImg + 'final_gjState' + '.png'
            plt.savefig(savename,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #---------Animations---------------------------------------------------------------------------------------

    if p.ani_ip32d is True and p.Ca_dyn is True and p.createAnimations is True:
        IP3plotting = np.asarray(sim.cIP3_time)
        IP3plotting = np.multiply(IP3plotting,1e3)

        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=IP3plotting,
            tit='IP3 concentration',
            cbtit='Concentration [umol/L]',
            clrAutoscale=p.autoscale_IP3_ani,
            clrMin=p.IP3_ani_min_clr,
            clrMax=p.IP3_ani_max_clr,
            save=p.saveAnimations,
            number_cells=p.enumerate_cells,
            saveFolder='animation/IP3',
            saveFile='ip3_',
            ignore_simECM=True,
            current_overlay=p.I_overlay,
        )

    if p.ani_dye2d is True and p.voltage_dye == 1 and \
       p.createAnimations is True:
        if p.sim_ECM is False:
            Dyeplotting = np.asarray(sim.cDye_time)
            Dyeplotting = np.multiply(Dyeplotting,1e3)

            AnimateCellData(
                sim=sim, cells=cells, p=p,
                zdata_t=Dyeplotting,
                tit='Morphogen Concentration',
                cbtit='Concentration [umol/L]',
                clrAutoscale=p.autoscale_Dye_ani,
                clrMin=p.Dye_ani_min_clr,
                clrMax=p.Dye_ani_max_clr,
                save=p.saveAnimations,
                number_cells=p.enumerate_cells,
                saveFolder='animation/Morphogen',
                saveFile='morphogen_',
                ignore_simECM=True,
            )

        else:
            AnimateDyeData(
                sim,cells,p,
                save=p.saveAnimations,
                ani_repeat=True,
                current_overlay=p.I_overlay,
                clrAutoscale=p.autoscale_Dye_ani,
                clrMin=p.Dye_ani_min_clr,
                clrMax=p.Dye_ani_max_clr,
                clrmap=p.default_cm,
                number_cells=p.enumerate_cells,
                saveFolder='animation/Dye',
                saveFile='Dye_',
            )

    if p.ani_ca2d is True and p.ions_dict['Ca'] == 1 and \
       p.createAnimations is True:
        tCa = [1e6*arr[sim.iCa] for arr in sim.cc_time]

        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=tCa,
            tit='Cytosolic Ca2+',
            cbtit='Concentration [nmol/L]',
            save=p.saveAnimations,
            clrAutoscale=p.autoscale_Ca_ani,
            clrMin=p.Ca_ani_min_clr,
            clrMax=p.Ca_ani_max_clr,
            number_cells=p.enumerate_cells,
            saveFolder='animation/Ca',
            saveFile='ca_',
            ignore_simECM=True,
        )

    if p.ani_pH2d is True and p.ions_dict['H'] == 1 and \
       p.createAnimations is True:
        tpH = [-np.log10(arr[sim.iH]) for arr in sim.cc_time]

        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=tpH,
            tit='Cytosolic pH',
            cbtit='pH',
            save=p.saveAnimations,
            clrAutoscale=p.autoscale_Ca_ani,
            clrMin=p.Ca_ani_min_clr,
            clrMax=p.Ca_ani_max_clr,
            number_cells=p.enumerate_cells,
            saveFolder='animation/pH',
            saveFile='pH_',
            ignore_simECM=True,
        )

    if p.ani_vm2d is True and p.createAnimations is True:
        if p.sim_ECM is False:
            vmplt = [1000*arr for arr in sim.vm_time]
        else:
            vmplt = [1000*arr for arr in sim.vm_Matrix]

        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=vmplt,
            tit='Cell Vmem',
            cbtit='Voltage [mV]',
            save=p.saveAnimations,
            clrAutoscale=p.autoscale_Vmem_ani,
            clrMin=p.Vmem_ani_min_clr,
            clrMax=p.Vmem_ani_max_clr,
            number_cells=p.enumerate_cells,
            current_overlay=p.I_overlay,
            saveFolder='animation/Vmem',
            saveFile='vm_',
            ignore_simECM=False,
        )

    if p.ani_vmgj2d is True and p.createAnimations is True:
        AnimateGJData(
            cells, sim, p,
            tit='Vcell ',
            save=p.saveAnimations,
            ani_repeat=True,
            saveFolder='animation/Vmem_gj',
            clrAutoscale=p.autoscale_Vgj_ani,
            clrMin=p.Vgj_ani_min_clr,
            clrMax=p.Vgj_ani_max_clr,
            saveFile='vmem_gj_',
            number_cells=False,
        )

    if p.ani_vcell is True and p.createAnimations is True and p.sim_ECM == 1:
        vcellplt = [1000*arr for arr in sim.vcell_time]

        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=vcellplt,
            tit='V in cell',
            cbtit='Voltage [mV]',
            clrAutoscale=p.autoscale_vcell_ani,
            clrMin=p.vcell_ani_min_clr,
            clrMax=p.vcell_ani_max_clr,
            save=p.saveAnimations,
            number_cells=p.enumerate_cells,
            saveFolder='animation/vcell',
            saveFile='vcell_',
            ignore_simECM=True,
            current_overlay=p.I_overlay,
        )

    if p.ani_I is True and p.createAnimations is True:
        AnimateCurrent(
            sim,cells,time,p,
            save=p.saveAnimations,
            ani_repeat=True,
            current_overlay=p.I_overlay,
            gj_current=True,
            clrAutoscale=p.autoscale_I_ani,
            clrMin=p.I_ani_min_clr,
            clrMax=p.I_ani_max_clr,
            clrmap=p.background_cm,
            number_cells=False,
            saveFolder='animation/current_gj',
            saveFile='I_',
        )

        if p.sim_ECM is True:
            AnimateCurrent(
                sim,cells,time,p,
                save=p.saveAnimations,
                ani_repeat=True,
                current_overlay=p.I_overlay,
                gj_current=False,
                clrAutoscale=p.autoscale_I_ani,
                clrMin=p.I_ani_min_clr,
                clrMax=p.I_ani_max_clr,
                clrmap=p.background_cm,
                number_cells=False,
                saveFolder='animation/current_ecm',
                saveFile='I_',
            )

    if p.ani_Efield is True and p.createAnimations is True:
        Ex_gj = sim.efield_gj_x_time
        Ey_gj = sim.efield_gj_y_time

        AnimateField(
            Ex_gj,Ey_gj,sim,cells,p,
            ani_repeat=True,
            save=p.saveAnimations,
            plot_ecm=False,
            title="Electric Field",
            cb_title="Electric Field [V/m]",
            colorAutoscale=p.autoscale_Efield_ani,
            colorMin=p.Efield_ani_min_clr,
            colorMax=p.Efield_ani_max_clr,
            saveFolder='animation/Efield',
            saveFile='Efield_',
        )

        if p.sim_ECM is True:

            Ex_ecm = sim.efield_ecm_x_time
            Ey_ecm = sim.efield_ecm_y_time

            AnimateField(
                Ex_ecm,Ey_ecm,sim,cells,p,
                ani_repeat=True,
                save=p.saveAnimations,
                plot_ecm=True,
                title="Electric Field",
                cb_title="Electric Field [V/m]",
                colorAutoscale=p.autoscale_Efield_ani,
                colorMin=p.Efield_ani_min_clr,
                colorMax=p.Efield_ani_max_clr,
                saveFolder='animation/Efield',
                saveFile='Efield_',
            )

    if p.ani_Velocity is True and p.fluid_flow is True and \
       p.deform_electro is True and p.createAnimations is True:
        AnimateVelocity(sim, cells, p, ani_repeat=True, save=p.saveAnimations)

    if p.ani_Deformation is True and p.deformation is True and \
       p.createAnimations is True and sim.run_sim is True:
        AnimateDeformation(sim, cells, p, ani_repeat=True, save=p.saveAnimations)

    if p.sim_eosmosis is True and p.sim_ECM is True and \
       cells.gradMem is not None:
        viz.plotMemData(cells,p,zdata=sim.rho_pump,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimension [um]')
        plt.title('Membrane ion pump density factor')

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        viz.plotMemData(cells,p,zdata=sim.rho_channel,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimension [um]')
        plt.title('Membrane ion channel density factor')

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    if p.ani_Pcell is True and p.deform_osmo is True and \
       p.createAnimations is True:
        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=sim.P_cells_time,
            tit='Hydrostatic Pressure in Cells',
            cbtit='Pressure [Pa]',
            clrAutoscale=p.autoscale_Pcell_ani,
            clrMin=p.Pcell_ani_min_clr,
            clrMax=p.Pcell_ani_max_clr,
            save=p.saveAnimations,
            number_cells=p.enumerate_cells,
            saveFolder='animation/Pcell',
            saveFile='Pcell_',
            ignore_simECM=True,
            current_overlay=p.I_overlay,
        )

    if p.ani_Pcell is True and p.createAnimations is True and \
       p.deform_osmo is True:
        AnimateCellData(
            sim=sim, cells=cells, p=p,
            zdata_t=sim.osmo_P_delta_time,
            tit='Osmotic Pressure in Cells',
            cbtit='Pressure [Pa]',
            clrAutoscale=p.autoscale_Pcell_ani,
            clrMin=p.Pcell_ani_min_clr,
            clrMax=p.Pcell_ani_max_clr,
            save=p.saveAnimations,
            number_cells=p.enumerate_cells,
            saveFolder='animation/OsmoP',
            saveFile='OsmoP_',
            ignore_simECM=True,
            current_overlay=p.I_overlay,
        )

    if p.ani_force is True and p.createAnimations is True:
        if p.deform_electro is True:
            FEx = [(1/p.um)*arr for arr in sim.F_electro_x_time]
            FEy = [(1/p.um)*arr for arr in sim.F_electro_y_time]

            AnimateField(
                FEx,FEy,sim,cells,p,
                ani_repeat=True,
                save=p.saveAnimations,
                plot_ecm=False,
                saveFolder='animation/ElectrostaticFfield',
                saveFile='EFfield_',
                title="Electrostatic Body Force",
                cb_title="Force [N/cm3]",
                colorAutoscale=p.autoscale_force_ani,
                colorMin=p.force_ani_min_clr,
                colorMax=p.force_ani_max_clr,
            )

        if p.deform_osmo is True:
            FHx = [(1/p.um)*arr for arr in sim.F_hydro_x_time]
            FHy = [(1/p.um)*arr for arr in sim.F_hydro_y_time]

            AnimateField(
                FHx,FHy,sim,cells,p,
                ani_repeat=True,
                save=p.saveAnimations,
                plot_ecm=False,
                saveFolder='animation/HydroFfield',
                saveFile='OFfield_',
                title="Hydrostatic Body Force",
                cb_title="Force [N/cm3]",
                colorAutoscale=p.autoscale_force_ani,
                colorMin=p.force_ani_min_clr,
                colorMax=p.force_ani_max_clr,
            )

    if p.ani_venv is True and p.createAnimations is True and p.sim_ECM is True:
        AnimateEnv(
            sim,cells,sim.time,p,
            clrAutoscale=p.autoscale_venv_ani,
            clrMin=p.venv_min_clr,
            clrMax=p.venv_max_clr,
            save=p.saveAnimations,
        )

    if p.ani_mem is True and p.sim_eosmosis is True and p.sim_ECM is True:
        AnimateMem(
            sim,cells,sim.time,p,
            clrAutoscale=p.autoscale_mem_ani,
            clrMin=p.mem_ani_min_clr,
            clrMax=p.mem_ani_max_clr,
            save=p.saveAnimations,
            current_overlay=p.I_overlay,
        )

    if p.turn_all_plots_off is False:
        plt.show()

    else:
        loggers.log_info(
            'As the config file results option "turn all plots off" is set to "True",\n'
            'plots and data have been exported to the results folder defined in the config\n'
            'file. To show animations, consider setting "turn all plots off" to "False".'
        )
