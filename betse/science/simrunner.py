#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


import os
import os.path
import time

from scipy import interpolate as interp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from betse.science import visualize as viz
from betse.science import filehandling as fh
from betse.science.sim import Simulator
from betse.science.parameters import Parameters
from betse.science.cells import Cells
from betse.science.tissue.handler import TissueHandler
from betse.util.io import loggers
from betse.util.path import files, paths
from betse.exceptions import BetseExceptionSimulation, BetseExceptionParameters


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

            if p.fluid_flow is True: # if user desires fluid flow:

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
            if p.fluid_flow is True:
                loggers.log_info('Creating cell network Poisson solver...')
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

            elif p.autoInit is False:
                raise BetseExceptionSimulation("Run terminated due to missing seed.\n "
                                               "Please run 'betse seed' to try again.")

        sim = Simulator(p)   # create an instance of Simulator

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
            loggers.log_info('When ready, close all of the figure windows to proceed with scheduled simulation runs.')
            plots4Sim(
                p.plot_cell,cells,sim,p,
                saveImages = p.autosave,
                animate=p.createAnimations,
                saveAni=p.saveAnimations)
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

        if p.sim_ECM is False:

            sim.runSim(cells,p,save=True)   # run and optionally save the simulation to the cache

        elif p.sim_ECM is True:

            sim.runSim_ECM(cells,p,save=True)   # run and optionally save the simulation to the cache

        loggers.log_info(
            'The simulation took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        loggers.log_info('When ready, close all of the figure windows to end the program.')

        if p.turn_all_plots_off is False:
            plots4Sim(
                p.plot_cell,cells,sim,p,
                saveImages = p.autosave,
                animate=p.createAnimations,
                saveAni=p.saveAnimations)
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
            raise BetseExceptionSimulation("Ooops! No such initialization file found to plot!")

        plots4Sim(p.plot_cell,cells,sim,p,saveImages=p.autosave, animate=p.createAnimations,
            saveAni=p.saveAnimations)
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
            raise BetseExceptionSimulation("Ooops! No such simulation file found to plot!")

        plots4Sim(
            p.plot_cell,cells,sim,p,
            saveImages=p.autosave,
            animate=p.createAnimations,
            saveAni=p.saveAnimations)
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

            images_path = p.sim_results
            image_cache_dir = os.path.expanduser(images_path)
            os.makedirs(image_cache_dir, exist_ok=True)
            savedImg = os.path.join(image_cache_dir, 'fig_')


        fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(p,dyna,cells)

        if p.autosave is True:
            savename10 = savedImg + 'cluster_mosaic' + '.png'
            plt.savefig(savename10,format='png')

        plt.show(block = False)

        if p.sim_ECM is True:

            plt.figure()
            ax99 = plt.subplot(111)
            plt.imshow(np.log10(sim.D_env_weight.reshape(cells.X.shape)),origin='lower',
                extent= [p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],cmap=p.default_cm)
            plt.colorbar()

            cell_edges_flat = cells.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(1.0)
            ax99.add_collection(coll)

            plt.title('Logarithm of Environmental Diffusion Weight Matrix')

            if p.autosave is True:
                savename10 = savedImg + 'env_diffusion_weights' + '.png'
                plt.savefig(savename10,format='png')

            plt.show(block =False)

            plt.figure()
            plt.imshow(cells.maskM,origin='lower',
                       extent= [p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax])
            plt.colorbar()
            plt.title('Cluster Masking Matrix')
            plt.show(block=False)

        # plot gj
        fig_x = plt.figure()
        ax_x = plt.subplot(111)

        # cell_edges_flat = p.um*cells.mem_edges_flat
        # coll = LineCollection(cell_edges_flat,color='k',linewidth=1.0)
        # ax_x.add_collection(coll)


        con_segs = cells.nn_edges
        connects = p.um*np.asarray(con_segs)
        collection = LineCollection(connects,linewidths=2.0,color='b')
        ax_x.add_collection(collection)
        plt.axis('equal')
        plt.axis([cells.xmin*p.um,cells.xmax*p.um,cells.ymin*p.um,cells.ymax*p.um])

        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Cell Connectivity Network')

        if p.autosave is True:
            savename10 = savedImg + 'gj_connectivity_network' + '.png'
            plt.savefig(savename10,format='png')

        plt.show()


def plots4Sim(plot_cell,cells,sim,p, saveImages=False, animate=0,saveAni=False):
    if saveImages is True:
        images_path = p.sim_results
        image_cache_dir = os.path.expanduser(images_path)
        os.makedirs(image_cache_dir, exist_ok=True)
        savedImg = os.path.join(image_cache_dir, 'fig_')

    #--------Single cell data graphs-----------------------------------------------------------------------------------

    if p.plot_single_cell_graphs is True:
        figConcsNa, axConcsNa = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,plot_cell,fig=None,
             ax=None,lncolor='g',ionname='Na+')

        titNa = 'Sodium concentration in cell ' + str(plot_cell)
        axConcsNa.set_title(titNa)

        if saveImages is True:
            savename1 = savedImg + 'concNa_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figConcsK, axConcsK = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,plot_cell,fig=None,
            ax=None,lncolor='b',ionname='K+')

        titK = 'Potassium concentration in cell ' + str(plot_cell)
        axConcsK.set_title(titK)

        if saveImages is True:
            savename1 = savedImg + 'concK_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figConcsM, axConcsM = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,plot_cell,fig=None,
             ax=None,lncolor='r',ionname='M-')

        titM = 'M Anion concentration in cell ' + str(plot_cell)
        axConcsM.set_title(titM)

        if saveImages is True:
            savename1 = savedImg + 'concM_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim,plot_cell,p,fig=None,ax=None,lncolor='k')
        titV = 'Voltage (Vmem) in cell ' + str(plot_cell)
        axVt.set_title(titV)

        if saveImages is True:
            savename2 = savedImg + 'Vmem_time' + '.png'
            plt.savefig(savename2,dpi=300,format='png')

        plt.show(block=False)

        #-----plot single cell transmembrane current-------------
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
        axI.set_title('Transmembrane current density for cell ' + str(plot_cell) )
        axI.set_xlabel('Time [s]')
        axI.set_ylabel('Current density [A/m2]')

        if saveImages is True:
            savename = savedImg + 'Imem_time' + '.png'
            plt.savefig(savename,dpi=300,format='png')

        plt.show(block=False)


        #---------------------------------------------------------

        if p.ions_dict['Ca'] ==1:
            figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,plot_cell,fig=None,
                 ax=None,lncolor='g',ionname='Ca2+ cell')
            titCa =  'Cytosolic Ca2+ in cell index ' + str(plot_cell)
            axA.set_title(titCa)

            if saveImages is True:
                savename3 = savedImg + 'cytosol_Ca_time' + '.png'
                plt.savefig(savename3,dpi=300,format='png')

            plt.show(block=False)

            if p.Ca_dyn == 1:
                figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,plot_cell,fig=None,
                     ax=None,lncolor='b',ionname='Ca2+ cell')
                titER =  'ER Ca2+ in cell index ' + str(plot_cell)
                axD.set_title(titER)

                if saveImages is True:
                    savename4 = savedImg + 'ER_Ca_time' + '.png'
                    plt.savefig(savename4,dpi=300,format='png')

                plt.show(block=False)

                figPro, axPro = viz.plotSingleCellData(sim.time, sim.cIP3_time,plot_cell, lab='IP3 [mmol/L]')
                titIP3 =  'IP3 in cell index ' + str(plot_cell)
                axPro.set_title(titIP3)

                if saveImages is True:
                    savename5 = savedImg + 'IP3_time' + '.png'
                    plt.savefig(savename5,dpi=300,format='png')

                plt.show(block=False)

    #-----2D data graphs-----------------------------------------------------------------------------------------------

    if p.plot_venv is True and p.sim_ECM is True:

        plt.figure()
        venv_plt = plt.imshow(1000*sim.v_env.reshape(cells.X.shape),origin='lower',
            extent= [p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],cmap=p.default_cm)
        if p.autoscale_venv is False:

            venv_plt.set_clim(p.venv_min_clr, p.venv_max_clr)

        plt.colorbar()
        plt.title('Environmental Voltage [mV]')

        if saveImages is True:
            savename10 = savedImg + 'Final_environmental_V' + '.png'
            plt.savefig(savename10,format='png')

        plt.show(block=False)

    if p.plot_rho2d is True:

        if p.data_type_rho == 'ECM' and p.sim_ECM is True:

            plt.figure()
            plt.imshow(sim.rho_env.reshape(cells.X.shape)*(1/p.ff_env),origin='lower',
                extent= [p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],cmap=p.default_cm)
            plt.colorbar()
            plt.title('Environmental Charge Density [C/m3]')

            if saveImages is True:
                savename10 = savedImg + 'Final_environmental_charge' + '.png'
                plt.savefig(savename10,format='png')

            plt.show(block =False)

        elif p.data_type_rho == 'GJ':

            if p.showCells is True:


                figX, axX, cbX = viz.plotPolyData(sim,cells,p,zdata=(sim.rho_cells)*(1/p.ff_cell),number_cells=p.enumerate_cells,
                    clrAutoscale = p.autoscale_rho, clrMin = p.rho_min_clr, clrMax = p.rho_max_clr,
                    clrmap = p.default_cm,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)

            else:

                figX, axX, cbX = viz.plotCellData(sim,cells,p,zdata = sim.rho_cells*(1/p.ff_cell),clrAutoscale = p.autoscale_rho,
                        clrMin = p.rho_min_clr, clrMax = p.rho_max_clr, clrmap = p.default_cm,
                        number_cells=p.enumerate_cells, current_overlay=p.I_overlay,plotIecm=p.IecmPlot)

            figX.suptitle('Final Cell Charge Density',fontsize=14, fontweight='bold')
            axX.set_xlabel('Spatial distance [um]')
            axX.set_ylabel('Spatial distance [um]')
            cbX.set_label('Net Charge Density [C/m3]')

            if saveImages is True:
                savename9 = savedImg + 'final_cellCharge' + '.png'
                plt.savefig(savename9,format='png')

            plt.show(block=False)

    if p.plot_vcell2d is True and p.sim_ECM is True:

        if p.showCells is True:

            figX, axX, cbX = viz.plotPolyData(sim,cells,p,zdata=sim.vcell_time[-1]*1e3,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_vcell, clrMin = p.vcell_min_clr, clrMax = p.vcell_max_clr,
                clrmap = p.default_cm,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)

        else:

            figX, axX, cbX = viz.plotCellData(sim,cells,p,zdata=1000*sim.vcell_time[-1],clrAutoscale = p.autoscale_vcell,
                    clrMin = p.vcell_min_clr, clrMax = p.vcell_max_clr, clrmap = p.default_cm,
                    number_cells=p.enumerate_cells, current_overlay=p.I_overlay,plotIecm=p.IecmPlot)

        figX.suptitle('Final Cell Voltage',fontsize=14, fontweight='bold')
        axX.set_xlabel('Spatial distance [um]')
        axX.set_ylabel('Spatial distance [um]')
        cbX.set_label('Voltage mV')

        if saveImages is True:
            savename9 = savedImg + 'final_cellVoltage' + '.png'
            plt.savefig(savename9,format='png')

        plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_vm2d is True:

        if p.sim_ECM is True:

            figV, axV, cbV = viz.plotHetMem(sim,cells,p,zdata=1000*sim.vm_Matrix[-1],number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm,
                edgeOverlay = p.showCells,  number_ecm = p.enumerate_cells,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)


        elif p.sim_ECM is False:

            if p.showCells is True:
                figV, axV, cbV = viz.plotPolyData(sim,cells,p,zdata=1000*sim.vm_time[-1],clrAutoscale = p.autoscale_Vmem,
                    clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, number_cells=p.enumerate_cells,
                    clrmap = p.default_cm,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)
            else:

                figV, axV, cbV = viz.plotCellData(sim,cells,p,zdata=1000*sim.vm_time[-1],clrAutoscale = p.autoscale_Vmem,
                    clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm,
                    number_cells=p.enumerate_cells, current_overlay=p.I_overlay,plotIecm=p.IecmPlot)


        figV.suptitle('Final Vmem',fontsize=14, fontweight='bold')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        if saveImages is True:
            savename5 = savedImg + 'final_Vmem_2D' + '.png'
            plt.savefig(savename5,format='png')

        plt.show(block=False)

    if p.GHK_calc is True:

        if p.showCells is True:
            figV_ghk, axV_ghk, cbV_ghk = viz.plotPolyData(sim,cells,p,zdata=1000*sim.vm_GHK_time[-1],
                clrAutoscale = p.autoscale_Vmem,
                clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, number_cells=p.enumerate_cells,
                clrmap = p.default_cm,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)
        else:

            figV_ghk, axV_ghk, cbV_ghk = viz.plotCellData(sim,cells,p,zdata=1000*sim.vm_GHK_time[-1],
                clrAutoscale = p.autoscale_Vmem,
                clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm,
                number_cells=p.enumerate_cells, current_overlay=p.I_overlay,plotIecm=p.IecmPlot)

        figV_ghk.suptitle('Final Vmem using Goldman Equation',fontsize=14, fontweight='bold')
        axV_ghk.set_xlabel('Spatial distance [um]')
        axV_ghk.set_ylabel('Spatial distance [um]')
        cbV_ghk.set_label('Voltage [mV]')

        if saveImages is True:
            savename5 = savedImg + 'final_Vmem_GHK_2D' + '.png'
            plt.savefig(savename5,format='png')

        plt.show(block=False)


    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_ip32d is True and p.scheduled_options['IP3'] != 0 or p.Ca_dyn != 0:

        if p.sim_ECM is False:

            if p.showCells is True:
                figIP3, axIP3, cbIP3 = viz.plotPolyData(sim, cells,p,zdata=sim.cIP3_time[-1]*1e3,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_IP3, clrMin = p.IP3_min_clr, clrMax = p.IP3_max_clr, clrmap = p.default_cm)
            else:
                 figIP3, axIP3, cbIP3 = viz.plotCellData(sim,cells,p,zdata=sim.cIP3_time[-1]*1e3,number_cells=p.enumerate_cells,
                 clrAutoscale = p.autoscale_IP3, clrMin = p.IP3_min_clr, clrMax = p.IP3_max_clr, clrmap = p.default_cm)

        else:

            figIP3 = plt.figure()
            axIP3 = plt.subplot(111)

            ip3Env = sim.cIP3_env*1e3
            ip3Cell = sim.cIP3*1e3

            bkgPlot = axIP3.imshow(ip3Env.reshape(cells.X.shape),origin='lower',extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,
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

        if saveImages is True:
            savename6 = savedImg + 'final_IP3_2D' + '.png'
            plt.savefig(savename6,format='png')

        plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_dye2d is True and p.voltage_dye == 1:
        if p.sim_ECM is False:
            if p.showCells is True:
                figVdye, axVdye, cbVdye = viz.plotPolyData(sim, cells,p,zdata=sim.cDye_time[-1]*1e3,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_Dye, clrMin = p.Dye_min_clr, clrMax = p.Dye_max_clr, clrmap = p.default_cm)

            else:
                figVdye, axVdye, cbVdye = viz.plotCellData(sim,cells,p,zdata=sim.cDye_time[-1]*1e3,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_Dye, clrMin = p.Dye_min_clr, clrMax = p.Dye_max_clr, clrmap = p.default_cm)

        else:
            figVdye = plt.figure()
            axVdye = plt.subplot(111)

            dyeEnv = sim.cDye_env*1e3
            dyeCell = sim.cDye_cell*1e3

            bkgPlot = axVdye.imshow(dyeEnv.reshape(cells.X.shape),origin='lower',extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,
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
                cbVDye = figVdye.colorbar(coll)

            xmin = cells.xmin*p.um
            xmax = cells.xmax*p.um
            ymin = cells.ymin*p.um
            ymax = cells.ymax*p.um

            axVdye.axis([xmin,xmax,ymin,ymax])

        axVdye.set_title('Final Morphogen Concentration')
        axVdye.set_xlabel('Spatial distance [um]')
        axVdye.set_ylabel('Spatial distance [um]')
        # cbVdye.set_label('Concentration umol/L')

        if saveImages is True:
            savename7 = savedImg + 'final_morphogen_2D' + '.png'
            plt.savefig(savename7,format='png')

        plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_ca2d is True and p.ions_dict['Ca'] == 1:

        if p.showCells is True:
            figCa, axCa, cbCa = viz.plotPolyData(sim,cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,number_cells= p.enumerate_cells,
            clrAutoscale = p.autoscale_Ca, clrMin = p.Ca_min_clr, clrMax = p.Ca_max_clr, clrmap = p.default_cm)
        else:
            figCa, axCa, cbCa = viz.plotCellData(sim,cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,number_cells=p.enumerate_cells,
            clrAutoscale = p.autoscale_Ca, clrMin = p.Ca_min_clr, clrMax = p.Ca_max_clr, clrmap = p.default_cm)

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if saveImages is True:
            savename8 = savedImg + 'final_Ca_2D' + '.png'
            plt.savefig(savename8,format='png')

        plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_I2d is True:

        figI, axI, cbI = viz.streamingCurrent(sim, cells,p,plot_Iecm = False, clrAutoscale = p.autoscale_I2d,
            clrMin = p.I_min_clr, clrMax = p.I_max_clr, clrmap= p.default_cm,
            edgeOverlay = p.showCells,number_cells = p.enumerate_cells)

        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current Density [A/m2]')

        if saveImages is True:
            savename10 = savedImg + 'Final_Current_gj' + '.png'
            plt.savefig(savename10,format='png')

        plt.show(block=False)

        if p.sim_ECM is True:

            figI2, axI2, cbI2 = viz.streamingCurrent(sim, cells,p,plot_Iecm = True, clrAutoscale = p.autoscale_I2d,
            clrMin = p.I_min_clr, clrMax = p.I_max_clr, clrmap= p.default_cm,
            edgeOverlay = p.showCells,number_cells = p.enumerate_cells)

            axI2.set_xlabel('Spatial distance [um]')
            axI2.set_ylabel('Spatial distance [um]')
            cbI2.set_label('Current Density [A/m2]')

            if saveImages is True:
                savename11 = savedImg + 'Final_Current_extracellular' + '.png'
                plt.savefig(savename11,format='png')

            plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_Efield is True:

        viz.plotEfield(sim,cells,p)

        if saveImages is True:
            savename12 = savedImg + 'Final_Electric_Field' + '.png'
            plt.savefig(savename12,format='png')

        plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_P is True:

        if p.fluid_flow is True or p.deformation is True:

            if p.showCells is True:
                figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=sim.P_cells,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_P, clrMin = p.P_min_clr, clrMax = p.P_max_clr, clrmap = p.default_cm)
            else:
                 figP, axP, cbP = viz.plotCellData(sim,cells,p,zdata=sim.P_cells,number_cells=p.enumerate_cells,
                 clrAutoscale = p.autoscale_P, clrMin = p.P_min_clr, clrMax = p.P_max_clr, clrmap = p.default_cm)

            axP.set_title('Final Hydrostatic Pressure in Cell Network')
            axP.set_xlabel('Spatial distance [um]')
            axP.set_ylabel('Spatial distance [um]')
            cbP.set_label('Pressure [Pa]')

            if saveImages is True:
                savename13 = savedImg + 'final_P_2D_gj' + '.png'
                plt.savefig(savename13,format='png')

            plt.show(block=False)


            if p.sim_ECM is True and p.fluid_flow is True:

                plt.figure()
                plt.imshow(sim.P_env,origin='lower',extent=[cells.xmin,cells.xmax,cells.ymin,cells.ymax],cmap=p.default_cm)
                plt.colorbar()
                plt.axis('equal')
                plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
                plt.title('Final Extracellular Pressure [Pa]')

                if saveImages is True:
                    savename13 = savedImg + 'final_P_2D_env' + '.png'
                    plt.savefig(savename13,format='png')

                plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_osmoP is True: # FIXME make this into a pressure head plot set

        if p.deform_osmo is True:

            if p.showCells is True:

                osmo_P = sim.osmo_P_delta

                # P_cell = np.dot(cells.M_sum_mems,sim.P_mem)/cells.num_mems

                figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=osmo_P,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr, clrmap = p.default_cm)

            else:
                 figP, axP, cbP = viz.plotCellData(sim,cells,p,zdata=osmo_P,number_cells=p.enumerate_cells,
                 clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr, clrmap = p.default_cm)


            axP.set_title('Final Osmotic Pressure in Cell Network')
            axP.set_xlabel('Spatial distance [um]')
            axP.set_ylabel('Spatial distance [um]')
            cbP.set_label('Pressure Difference Cell Interior vs Exterior [Pa]')


            if saveImages is True:
                savename13 = savedImg + 'final_osmoP_2D' + '.png'
                plt.savefig(savename13,format='png')

            plt.show(block=False)

        if p.gravity is True:

            # gravity hydrostatic pressure head:
            if p.showCells is True:

                gravP = np.dot(cells.M_sum_mems,sim.P_gravity)/cells.num_mems

                figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=gravP,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr, clrmap = p.default_cm)

            else:
                 figP, axP, cbP = viz.plotCellData(sim,cells,p,zdata=gravP,number_cells=p.enumerate_cells,
                 clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr, clrmap = p.default_cm)


            axP.set_title('Final Gravity-Induced Pressure Head in Cell Network')
            axP.set_xlabel('Spatial distance [um]')
            axP.set_ylabel('Spatial distance [um]')
            cbP.set_label('Pressure [Pa]')


            if saveImages is True:
                savename13 = savedImg + 'final_gravityP_2D' + '.png'
                plt.savefig(savename13,format='png')

            plt.show(block=False)

        #-----Electrostatic stress------------------------------------------------------------

        if p.deform_electro is True:

            # average the electrostatic pressure from mems to whole cell:
            P_electro_cell = np.dot(cells.M_sum_mems,sim.P_electro)/cells.num_mems

            figEP, axEP, cbEP = viz.plotPolyData(sim, cells,p,zdata=P_electro_cell,number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_osmoP, clrMin = p.osmoP_min_clr, clrMax = p.osmoP_max_clr, clrmap = p.default_cm)

            # axEP.quiver(cells.mem_vects_flat[:,0]*p.um,cells.mem_vects_flat[:,1]*p.um,
            #     sim.P_electro*cells.mem_vects_flat[:,2],sim.P_electro*cells.mem_vects_flat[:,3])

            axEP.set_title('Final Electrostatic Pressure in Cell Network')
            axEP.set_xlabel('Spatial distance [um]')
            axEP.set_ylabel('Spatial distance [um]')
            cbEP.set_label('Pressure [Pa]')


            if saveImages is True:
                savename13 = savedImg + 'final_electroP_2D' + '.png'
                plt.savefig(savename13,format='png')

            plt.show(block=False)


    if p.plot_Vel is True:

        if p.fluid_flow is True:

            ucellso = sim.u_cells_x
            vcellso = sim.u_cells_y

            ucells = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),
                                     ucellso,(cells.Xgrid,cells.Ygrid), fill_value=0)

            ucells = ucells*cells.maskM

            vcells = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),
                                     vcellso,(cells.Xgrid,cells.Ygrid), fill_value=0)

            vcells = vcells*cells.maskM

            Ucells = np.sqrt(ucells**2 + vcells**2)*1e9

            if Ucells.max() != 0.0:
                lw = (Ucells/Ucells.max()) + 0.5

            plt.figure()
            plt.imshow(Ucells,origin='lower',extent=[cells.xmin,cells.xmax,cells.ymin,cells.ymax],cmap=p.default_cm)
            plt.colorbar()
            plt.streamplot(cells.Xgrid,cells.Ygrid,ucells/Ucells.max(),vcells/Ucells.max(),density=p.stream_density,linewidth=lw,color='k')
            plt.axis('equal')
            plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
            plt.title('Final Fluid Velocity in Cell Collective [nm/s]')

            if saveImages is True:
                savename13 = savedImg + 'final_vel_2D_gj' + '.png'
                plt.savefig(savename13,format='png')

            plt.show(block=False)

            if p.sim_ECM is True and p.fluid_flow is True:

                u = sim.u_at_c
                v = sim.v_at_c
                U = np.sqrt(u**2 + v**2)*1e9

                plt.figure()
                plt.imshow(U,origin='lower',extent=[cells.xmin,cells.xmax,cells.ymin,cells.ymax],cmap=p.default_cm)
                plt.colorbar()
                plt.streamplot(cells.X,cells.Y,u,v,density=p.stream_density,color='k')
                plt.axis('equal')
                plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
                plt.title('Final Extracellular Fluid Velocity [nm/s]')

                if saveImages is True:
                    savename13 = savedImg + 'final_vel_2D_env' + '.png'
                    plt.savefig(savename13,format='png')

                plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.ani_ip32d is True and p.Ca_dyn is True and animate == 1:
        IP3plotting = np.asarray(sim.cIP3_time)
        IP3plotting = np.multiply(IP3plotting,1e3)

        if p.showCells is True:

            viz.AnimateCellData(sim,cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_IP3_ani, clrMin = p.IP3_ani_min_clr, clrMax = p.IP3_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/IP3',
                saveFile = 'ip3_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_IP3_ani, clrMin = p.IP3_ani_min_clr, clrMax = p.IP3_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/IP3', saveFile = 'ip3_')

    if p.ani_dye2d is True and p.voltage_dye == 1 and animate ==1:

        if p.sim_ECM is False:

            Dyeplotting = np.asarray(sim.cDye_time)
            Dyeplotting = np.multiply(Dyeplotting,1e3)


            if p.showCells is True:
                viz.AnimateCellData(sim,cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                    clrAutoscale = p.autoscale_Dye_ani, clrMin = p.Dye_ani_min_clr, clrMax = p.Dye_ani_max_clr, clrmap = p.default_cm,
                    save=saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Dye',
                    saveFile = 'dye_',ignore_simECM =True)
            else:
                viz.AnimateCellData_smoothed(sim,cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                    clrAutoscale = p.autoscale_Dye_ani, clrMin = p.Dye_ani_min_clr, clrMax = p.Dye_ani_max_clr, clrmap = p.default_cm,
                    save=saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/Dye', saveFile = 'Dye_')

        else:

            zenv_t = sim.cDye_env_time[:]

            viz.AnimateDyeData(sim,cells,p,save=saveAni,ani_repeat=True, current_overlay = p.I_overlay,
            clrAutoscale = p.autoscale_Dye_ani, clrMin = p.Dye_ani_min_clr, clrMax = p.Dye_ani_max_clr, clrmap = p.default_cm,
            number_cells = p.enumerate_cells, saveFolder = '/animation/Dye', saveFile = 'Dye_')

    if p.ani_ca2d is True and p.ions_dict['Ca'] == 1 and animate == 1:

        tCa = [1e6*arr[sim.iCa] for arr in sim.cc_time]

        if p.showCells is True:
            viz.AnimateCellData(sim,cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                clrAutoscale = p.autoscale_Ca_ani, clrMin = p.Ca_ani_min_clr, clrMax = p.Ca_ani_max_clr, clrmap = p.default_cm,
                ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Ca',
                saveFile = 'ca_',ignore_simECM = True)
        else:
            viz.AnimateCellData_smoothed(sim,cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                clrAutoscale = p.autoscale_Ca_ani, clrMin = p.Ca_ani_min_clr, clrMax = p.Ca_ani_max_clr, clrmap = p.default_cm,
                ani_repeat=True,number_cells=False,saveFolder = '/animation/Ca', saveFile = 'ca_')

    if p.ani_vm2d is True and animate == 1:

        vmplt = [1000*arr for arr in sim.vm_time]

        if p.sim_ECM is True:

            viz.AnimateCellData(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr,
                clrmap = p.default_cm, ani_repeat=True,number_cells=p.enumerate_cells, current_overlay=p.I_overlay,
                saveFolder = '/animation/Vmem', saveFile = 'vm_')

        elif p.sim_ECM is False:

            if p.showCells is True:
                viz.AnimateCellData(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                     clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr, clrmap = p.default_cm,
                    ani_repeat=True,number_cells=p.enumerate_cells, current_overlay=p.I_overlay,
                    saveFolder = '/animation/Vmem', saveFile = 'vm_')
            else:
                viz.AnimateCellData_smoothed(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                     clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr, clrmap = p.default_cm,
                    ani_repeat=True,number_cells=False,saveFolder = '/animation/Vmem', saveFile = 'vm_',current_overlay=p.I_overlay)

    if p.ani_vmgj2d is True and animate == 1:

        if p.sim_ECM is True:
            viz.AnimateGJData(cells, sim, p, tit='Vcell ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj_', number_cells=False)

        elif p.sim_ECM is False:

            if p.showCells is True:
                viz.AnimateGJData(cells, sim, p, tit='Vmem ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj_', number_cells=False)
            else:
                viz.AnimateGJData_smoothed(cells, sim, p, tit='Vmem ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj', number_cells=False)

    if p.ani_vcell is True and animate == 1 and p.sim_ECM == 1:

        vcellplt = [1000*arr for arr in sim.vcell_time]

        if p.showCells is True:

            viz.AnimateCellData(sim,cells,vcellplt,sim.time,p,tit='V in cell', cbtit = 'Voltage [mV]',
                clrAutoscale = p.autoscale_vcell_ani, clrMin = p.vcell_ani_min_clr, clrMax = p.vcell_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/vcell',
                saveFile = 'vcell_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,vcellplt,sim.time,p,tit='V in cell', cbtit = 'Voltage [mV]',
                clrAutoscale = p.autoscale_vcell_ani, clrMin = p.vcell_ani_min_clr, clrMax = p.vcell_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/vcell', saveFile = 'vcell_',
                current_overlay=p.I_overlay)

    if p.ani_I is True and animate == 1:

        viz.AnimateCurrent(sim,cells,time,p,save=saveAni,ani_repeat=True,current_overlay=p.I_overlay, gj_current =True,
            clrAutoscale=p.autoscale_I_ani,clrMin = p.I_ani_min_clr,clrMax = p.I_ani_max_clr,
            clrmap = p.default_cm, number_cells=False,saveFolder = '/animation/current_gj',saveFile = 'I_')

        if p.sim_ECM is True:

            viz.AnimateCurrent(sim,cells,time,p,save=saveAni,ani_repeat=True,current_overlay=p.I_overlay, gj_current =False,
            clrAutoscale=p.autoscale_I_ani,clrMin = p.I_ani_min_clr,clrMax = p.I_ani_max_clr,
            clrmap = p.default_cm, number_cells=False,saveFolder = '/animation/current_ecm',saveFile = 'I_')

    if p.ani_Efield is True and animate == 1:

        viz.AnimateEfield(sim,cells,p,ani_repeat = True, save = saveAni)

    if p.ani_Velocity is True and p.fluid_flow is True and animate == 1:

        viz.AnimateVelocity(sim,cells,p,ani_repeat = True, save = saveAni)

    if p.exportData is True:
        viz.exportData(cells, sim, p)


    if p.sim_eosmosis is True and p.sim_ECM is True and cells.gradMem is not None:

        viz.plotMemData(cells,p,zdata=sim.rho_pump,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimention [um]')
        plt.title('Membrane ion pump density factor')

        plt.show(block=False)

        viz.plotMemData(cells,p,zdata=sim.rho_channel,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimention [um]')
        plt.title('Membrane ion channel density factor')

        plt.show(block=False)

    if p.gj_flux_sensitive is True or p.v_sensitive_gj is True:

        # viz.plotMemData(cells,p,zdata=sim.rho_gj,clrmap=p.default_cm)
        fig_x = plt.figure()
        ax_x = plt.subplot(111)
        # plt.quiver(p.um*cells.nn_vects[:,0],p.um*cells.nn_vects[:,1],cells.nn_vects[:,2],cells.nn_vects[:,3],
        #     sim.gj_rho,cmap = p.default_cm)
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
        ax_x.set_title('Gap Junction Relative Permeability')

        plt.show(block=False)

    if p.ani_Pcell is True and animate == 1 and p.base_eosmo is True:

        if p.showCells is True:

            viz.AnimateCellData(sim,cells,sim.P_cells_time,sim.time,p,tit='Hydrostatic Pressure in Cells', cbtit = 'Pressure [Pa]',
                clrAutoscale = p.autoscale_Pcell_ani, clrMin = p.Pcell_ani_min_clr, clrMax = p.Pcell_ani_max_clr,
                clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Pcell',
                saveFile = 'Pcell_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,sim.P_cells_time,sim.time,p,tit='Hydrostatic Pressure in Cells', cbtit = 'Pressure [Pa]',
                clrAutoscale = p.autoscale_Pcell_ani, clrMin = p.Pcell_ani_min_clr, clrMax = p.Pcell_ani_max_clr,
                clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/Pcell', saveFile = 'Pcell_',
                current_overlay=p.I_overlay)

    if p.ani_osmoP is True and animate == 1:

        osmo_P_atm = [arr*(1) for arr in sim.osmo_P_delta_time]

        if p.showCells is True:

            viz.AnimateCellData(sim,cells,osmo_P_atm,sim.time,p,tit='Osmotic Pressure in Cells', cbtit = 'Pressure [Pa]',
                clrAutoscale = p.autoscale_osmoP_ani, clrMin = p.osmoP_ani_min_clr, clrMax = p.osmoP_ani_max_clr,
                clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/osmoP',
                saveFile = 'osmoP_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,osmo_P_atm,sim.time,p,tit='Osmotic Pressure in Cells',
                cbtit = 'Pressure [Pa]',
                clrAutoscale = p.autoscale_osmoP_ani, clrMin = p.osmoP_ani_min_clr, clrMax = p.osmoP_ani_max_clr,
                clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/osmoP', saveFile = 'osmoP_',
                current_overlay=p.I_overlay)

    if p.ani_venv is True and animate == 1 and p.sim_ECM is True:

        viz.AnimateEnv(sim,cells,sim.time,p,clrAutoscale=p.autoscale_venv_ani,clrMin=p.venv_min_clr,
                       clrMax=p.venv_max_clr, save = saveAni)

    if p.ani_mem is True and p.sim_eosmosis is True and p.sim_ECM is True:

        viz.AnimateMem(sim,cells,sim.time,p,clrAutoscale=p.autoscale_mem_ani,clrMin= p.mem_ani_min_clr,
                       clrMax=p.mem_ani_max_clr,save = saveAni,current_overlay=p.I_overlay)


    plt.show()


    #------------------------------------------------------------------------------------------------------------

    # Bx = sim.Bx
    # By = sim.By
    #
    # Bcx = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),Bx,(cells.X,cells.Y), fill_value=0)
    #
    # Bcy = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),By,(cells.X,cells.Y), fill_value=0)
    #
    # Bcells = np.sqrt(Bcx**2 + Bcy**2)
    #
    # plt.figure()
    # plt.imshow(Bcells,origin='lower',extent=[cells.xmin,cells.xmax,cells.ymin,cells.ymax],cmap=p.default_cm)
    # plt.colorbar()
    # plt.streamplot(cells.X,cells.Y,Bcx,Bcy,density=p.stream_density,color='k')
    # plt.axis('equal')
    # plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
    # plt.title('Final Magnetic Field in Cell Collective [unit]')

    # P = np.float64(sim.P_cells)
    #
    # plt.figure()
    # plt.tripcolor(cells.cell_centres[:,0],cells.cell_centres[:,1],P,shading='gouraud')
    # plt.colorbar()
    # plt.axis('equal')
    # plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
    # plt.title('Test Fluid Pressure in Cell Collective [um/s]')
    #
    # plt.show(block=False)
    #
    # ux = np.float64(sim.u_cells_x)
    # uy = np.float64(sim.u_cells_y)
    #
    # uxgj = np.float64(sim.u_cells_x_gj)
    # uygj = np.float64(sim.u_cells_y_gj)
    #
    # uu = np.sqrt(ux**2 + uy**2)*1e6
    #
    # plt.figure()
    # plt.tripcolor(cells.cell_centres[:,0],cells.cell_centres[:,1],uu,shading='gouraud')
    # plt.colorbar()
    # plt.quiver(cells.nn_vects[:,0],cells.nn_vects[:,1],uxgj,uygj)
    # plt.axis('equal')
    # plt.axis([cells.xmin,cells.xmax,cells.ymin,cells.ymax])
    # plt.title('Test Fluid Velocity in Cell Collective [um/s]')
    #
    # plt.show(block=False)


