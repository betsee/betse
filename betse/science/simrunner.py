#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


from betse.science import visualize as viz
from betse.science import filehandling as fh
from betse.science.compute import Simulator
from betse.science.parameters import Parameters
from betse.science.world import World
from betse.science.dynamics import Dynamics
from betse.util.io import loggers
from betse.util.path import files, paths
import matplotlib.pyplot as plt
import numpy as np
import os, os.path
import time
from betse.util.path import files
from betse.exceptions import BetseExceptionSimulation



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
        files.die_unless_found(config_filename)
        self._config_filename = config_filename
        self._config_basename = paths.get_basename(self._config_filename)

    def makeWorld(self,plotWorld = True):

        """
        In order to set up tissue profiles and other geometry-specific features, it is necessary
        to first create and plot the cells data structure. This will be loaded into the init
        and sim runs.

        """

        loggers.log_info(
            'Initializing simulation with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.I_overlay = False  # force the current overlay to be null
        sim = Simulator(p)   # create an instance of Simulator as it's needed by plotting objects

        if p.sim_ECM == False:

            cells = World(p,worldtype='basic')  # create an instance of world
            cells.containsECM = False
            loggers.log_info('Cell cluster is being created...')
            cells.makeWorld(p)     # call function to create the world

            # define the tissue and boundary profiles for plotting:
            loggers.log_info('Defining tissue and boundary profiles...')
            sim.baseInit(cells,p)
            dyna = Dynamics(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)

            cells.redo_gj(dyna,p)  # redo gap junctions to isolate different tissue types

            loggers.log_info('Cell cluster creation complete!')


            if plotWorld == True:

                fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(p,dyna,cells)

                plt.show(block=False)

                print(cells.cell_number)

        elif p.sim_ECM == True:

            cells = World(p,worldtype='full')  # create an instance of world
            cells.containsECM = True
            loggers.log_info('Cell cluster is being created...')
            cells.makeWorld(p)     # call function to create the world

            # define the tissue and boundary profiles for plotting:
            loggers.log_info('Defining tissue and boundary profiles...')
            sim.baseInit_ECM(cells,p)
            dyna = Dynamics(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)
            dyna.ecmBoundProfiles(sim,cells,p)

            cells.redo_gj(dyna,p)  # redo gap junctions to isolate different tissue types

            loggers.log_info('Cell cluster creation complete!')

            if plotWorld == True:

                fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(p,dyna,cells)

                plt.show(block=False)

                print(cells.cell_number)


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

        cells = World(p)  # create an instance of world


        if files.is_file(cells.savedWorld):
            cells,_ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            loggers.log_info('Cell cluster loaded.')

            if p.sim_ECM != cells.sim_ECM:
                loggers.log_info("Ooops! Cell cluster and config settings don't match!")
                loggers.log_info("Automatically creating cell cluster from current config file settings...")
                loggers.log_info("Warning: specified tissue profiles may no longer be correctly assigned.")
                self.makeWorld(plotWorld=False)  # create an instance of world
                loggers.log_info('Now using cell cluster to run initialization.')
                cells,_ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache


        else:
            loggers.log_info("Ooops! No such cell cluster file found to load!")

            if p.autoInit == True:
                loggers.log_info("Automatically creating cell cluster from config file settings...")
                self.makeWorld(plotWorld=False)  # create an instance of world
                loggers.log_info('Now using cell cluster to run initialization.')
                cells,_ = fh.loadWorld(cells.savedWorld)  # load the initialization from cache


            elif p.autoInit == False:
                raise BetseExceptionSimulation("Simulation terminated due to missing initialization. Please run"
                                               "an initialization and try again.")

        sim = Simulator(p)   # create an instance of Simulator

        if p.sim_ECM == False:

            sim.baseInit(cells, p)   # initialize simulation data structures
            # sim.tissueInit(cells,p)
            sim.runSim(cells,p)     # run and save the initialization

        elif p.sim_ECM == True:

            sim.baseInit_ECM(cells, p)   # initialize simulation data structures
            # sim.tissueInit(cells,p)
            sim.runSim_ECM(cells,p)     # run and save the initialization

        loggers.log_info('Initialization run complete!')
        loggers.log_info(
            'The initialization took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        if p.turn_all_plots_off == False:
            plots4Sim(
                p.plot_cell,cells,sim,p,
                saveImages = p.autosave,
                animate=p.createAnimations,
                saveAni=p.saveAnimations)
            plt.show()

        loggers.log_info('When ready, close all of the figure windows to proceed with scheduled simulation runs.')

        # if p.turn_all_plots_off == False:
        #     plots4Init(p.plot_cell,cells,sim,p,saveImages=p.autosave)
        #     plt.show()

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
            sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
            p.sim_ECM = cells.sim_ECM

        else:

            loggers.log_info("No initialization file found to run this simulation!")
            answer = p.autoInit

            if p.autoInit == True:
                loggers.log_info("Automatically running initialization...")
                self.initialize()
                loggers.log_info('Now using initialization to run simulation.')
                sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

            elif p.autoInit == False:
                raise BetseExceptionSimulation("Simulation terminated due to missing initialization. Please run"
                                               "an initialization and try again.")

        sim.fileInit(p)   # reinitialize save and load directories in case params defines new ones for this sim

        if p.sim_ECM == False:

            sim.runSim(cells,p,save=True)   # run and optionally save the simulation to the cache

        elif p.sim_ECM == True:

            sim.runSim_ECM(cells,p,save=True)   # run and optionally save the simulation to the cache

        loggers.log_info(
            'The simulation took {} seconds to complete.'.format(
                round(time.time() - start_time, 2)))

        loggers.log_info('When ready, close all of the figure windows to end the program.')

        if p.turn_all_plots_off == False:
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
            'Initializing simulation with configuration file "{}".'.format(
                self._config_basename))

        p = Parameters(config_filename = self._config_filename)     # create an instance of Parameters
        p.I_overlay = False # force the current overlay to be false as there's no data for it
        sim = Simulator(p)

        cells = World(p,worldtype='basic')

        if files.is_file(cells.savedWorld):
            cells,_ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
            p.sim_ECM = cells.sim_ECM
            loggers.log_info('Cell cluster loaded.')
        else:
            raise BetseExceptionSimulation("Ooops! No such cell cluster file found to load!")

        if p.sim_ECM == False:
            sim.baseInit(cells,p)
            dyna = Dynamics(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)
        else:
            sim.baseInit_ECM(cells,p)
            dyna = Dynamics(sim,cells,p)
            dyna.tissueProfiles(sim,cells,p)
            dyna.ecmBoundProfiles(sim,cells,p)



        fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(p,dyna,cells)

        plt.show()

def plots4Sim(plot_cell,cells,sim,p, saveImages=False, animate=0,saveAni=False):

    if saveImages == True:

        images_path = p.sim_results
        image_cache_dir = os.path.expanduser(images_path)
        os.makedirs(image_cache_dir, exist_ok=True)
        savedImg = os.path.join(image_cache_dir, 'fig_')

    if p.plot_single_cell_graphs == True:

        figConcsNa, axConcsNa = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,plot_cell,fig=None,
             ax=None,lncolor='g',ionname='Na+')

        titNa = 'Sodium concentration in cell ' + str(plot_cell)
        axConcsNa.set_title(titNa)

        if saveImages == True:
            savename1 = savedImg + 'concNa_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figConcsK, axConcsK = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,plot_cell,fig=None,
            ax=None,lncolor='b',ionname='K+')

        titK = 'Potassium concentration in cell ' + str(plot_cell)
        axConcsK.set_title(titK)

        if saveImages == True:
            savename1 = savedImg + 'concK_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figConcsM, axConcsM = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,plot_cell,fig=None,
             ax=None,lncolor='r',ionname='M-')

        titM = 'M Anion concentration in cell ' + str(plot_cell)
        axConcsM.set_title(titM)

        if saveImages == True:
            savename1 = savedImg + 'concM_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png')

        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,plot_cell,fig=None,ax=None,lncolor='b')
        titV = 'Voltage (Vmem) in cell ' + str(plot_cell)
        axVt.set_title(titV)

        if saveImages == True:
            savename2 = savedImg + 'Vmem_time' + '.png'
            plt.savefig(savename2,dpi=300,format='png')

        plt.show(block=False)

        if p.ions_dict['Ca'] ==1:
            figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,plot_cell,fig=None,
                 ax=None,lncolor='g',ionname='Ca2+ cell')
            titCa =  'Cytosolic Ca2+ in cell index ' + str(plot_cell)
            axA.set_title(titCa)

            if saveImages == True:
                savename3 = savedImg + 'cytosol_Ca_time' + '.png'
                plt.savefig(savename3,dpi=300,format='png')

            plt.show(block=False)

            if p.Ca_dyn == 1:
                figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,plot_cell,fig=None,
                     ax=None,lncolor='b',ionname='Ca2+ cell')
                titER =  'ER Ca2+ in cell index ' + str(plot_cell)
                axD.set_title(titER)

                if saveImages == True:
                    savename4 = savedImg + 'ER_Ca_time' + '.png'
                    plt.savefig(savename4,dpi=300,format='png')

                plt.show(block=False)

                figPro, axPro = viz.plotSingleCellData(sim.time, sim.cIP3_time,plot_cell, lab='IP3 [mmol/L]')
                titIP3 =  'IP3 in cell index ' + str(plot_cell)
                axPro.set_title(titIP3)

                if saveImages == True:
                    savename5 = savedImg + 'IP3_time' + '.png'
                    plt.savefig(savename5,dpi=300,format='png')

                plt.show(block=False)

    if p.plot_vcell2d == True and p.sim_ECM == True:

        if p.showCells == True:

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

        if saveImages == True:
            savename9 = savedImg + 'final_cellVoltage' + '.png'
            plt.savefig(savename9,format='png')

        plt.show(block=False)

    if p.plot_vm2d == True:

        if p.sim_ECM == True:

            figV, axV, cbV = viz.plotHetMem(sim,cells,p,zdata=1000*sim.vm_Matrix[-1],number_cells=p.enumerate_cells,
                clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm,
                edgeOverlay = p.showCells,  number_ecm = p.enumerate_cells,current_overlay = p.I_overlay,plotIecm=p.IecmPlot)

        elif p.sim_ECM == False:

            if p.showCells == True:
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

        if saveImages == True:
            savename5 = savedImg + 'final_Vmem_2D' + '.png'
            plt.savefig(savename5,format='png')

        plt.show(block=False)

    if p.plot_ip32d == True and p.scheduled_options['IP3'] != 0:

        if p.showCells == True:
            figIP3, axIP3, cbIP3 = viz.plotPolyData(sim, cells,p,zdata=sim.cIP3_time[-1]*1e3,number_cells=p.enumerate_cells,
            clrAutoscale = p.autoscale_IP3, clrMin = p.IP3_min_clr, clrMax = p.IP3_max_clr, clrmap = p.default_cm)
        else:
             figIP3, axIP3, cbIP3 = viz.plotCellData(sim,cells,p,zdata=sim.cIP3_time[-1]*1e3,number_cells=p.enumerate_cells,
             clrAutoscale = p.autoscale_IP3, clrMin = p.IP3_min_clr, clrMax = p.IP3_max_clr, clrmap = p.default_cm)

        axIP3.set_title('Final IP3 concentration')
        axIP3.set_xlabel('Spatial distance [um]')
        axIP3.set_ylabel('Spatial distance [um]')
        cbIP3.set_label('Concentration umol/L')

        if saveImages == True:
            savename6 = savedImg + 'final_IP3_2D' + '.png'
            plt.savefig(savename6,format='png')

        plt.show(block=False)

    if p.plot_dye2d == True and p.voltage_dye == 1:

        if p.showCells == True:
            figVdye, axVdye, cbVdye = viz.plotPolyData(sim, cells,p,zdata=sim.cDye_time[-1]*1e3,number_cells=p.enumerate_cells,
            clrAutoscale = p.autoscale_Dye, clrMin = p.Dye_min_clr, clrMax = p.Dye_max_clr, clrmap = p.default_cm)
        else:
            figVdye, axVdye, cbVdye = viz.plotCellData(sim,cells,p,zdata=sim.cDye_time[-1]*1e3,number_cells=p.enumerate_cells,
            clrAutoscale = p.autoscale_Dye, clrMin = p.Dye_min_clr, clrMax = p.Dye_max_clr, clrmap = p.default_cm)

        axVdye.set_title('Final voltage-sensitive dye')
        axVdye.set_xlabel('Spatial distance [um]')
        axVdye.set_ylabel('Spatial distance [um]')
        cbVdye.set_label('Concentration umol/L')

        if saveImages == True:
            savename7 = savedImg + 'final_dye_2D' + '.png'
            plt.savefig(savename7,format='png')

        plt.show(block=False)

    if p.plot_ca2d ==True and p.ions_dict['Ca'] == 1:

        if p.showCells == True:
            figCa, axCa, cbCa = viz.plotPolyData(sim,cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,number_cells= p.enumerate_cells,
            clrAutoscale = p.autoscale_Ca, clrMin = p.Ca_min_clr, clrMax = p.Ca_max_clr, clrmap = p.default_cm)
        else:
            figCa, axCa, cbCa = viz.plotCellData(sim,cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,number_cells=p.enumerate_cells,
            clrAutoscale = p.autoscale_Ca, clrMin = p.Ca_min_clr, clrMax = p.Ca_max_clr, clrmap = p.default_cm)

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if saveImages == True:
            savename8 = savedImg + 'final_Ca_2D' + '.png'
            plt.savefig(savename8,format='png')

        plt.show(block=False)

    if p.plot_I2d == True:

        figI, axI, cbI = viz.streamingCurrent(sim, cells,p,plot_Iecm = False, clrAutoscale = p.autoscale_I2d,
            clrMin = p.I_min_clr, clrMax = p.I_max_clr, clrmap= p.default_cm,
            edgeOverlay = p.showCells,number_cells = p.enumerate_cells)

        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current [pA]')

        if saveImages == True:
            savename10 = savedImg + 'Final_Current_gj' + '.png'
            plt.savefig(savename10,format='png')

        plt.show(block=False)

        if p.sim_ECM == True:

            figI2, axI2, cbI2 = viz.streamingCurrent(sim, cells,p,plot_Iecm = True, clrAutoscale = p.autoscale_I2d,
            clrMin = p.I_min_clr, clrMax = p.I_max_clr, clrmap= p.default_cm,
            edgeOverlay = p.showCells,number_cells = p.enumerate_cells)

            axI2.set_xlabel('Spatial distance [um]')
            axI2.set_ylabel('Spatial distance [um]')
            cbI2.set_label('Current [pA]')

            if saveImages == True:
                savename11 = savedImg + 'Final_Current_extracellular' + '.png'
                plt.savefig(savename11,format='png')

            plt.show(block=False)


    if p.ani_ip32d ==True and p.scheduled_options['IP3'] != 0 and animate == 1:
        IP3plotting = np.asarray(sim.cIP3_time)
        IP3plotting = np.multiply(IP3plotting,1e3)

        if p.showCells == True:

            viz.AnimateCellData(sim,cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_IP3_ani, clrMin = p.IP3_ani_min_clr, clrMax = p.IP3_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/IP3',
                saveFile = 'ip3_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_IP3_ani, clrMin = p.IP3_ani_min_clr, clrMax = p.IP3_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/IP3', saveFile = 'ip3_')

    if p.ani_dye2d == True and p.voltage_dye == 1 and animate ==1:

        Dyeplotting = np.asarray(sim.cDye_time)
        Dyeplotting = np.multiply(Dyeplotting,1e3)

        if p.showCells == True:
            viz.AnimateCellData(sim,cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_Dye_ani, clrMin = p.Dye_ani_min_clr, clrMax = p.Dye_ani_max_clr, clrmap = p.default_cm,
                save=saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Dye',
                saveFile = 'dye_',ignore_simECM =True)
        else:
            viz.AnimateCellData_smoothed(sim,cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                clrAutoscale = p.autoscale_Dye_ani, clrMin = p.Dye_ani_min_clr, clrMax = p.Dye_ani_max_clr, clrmap = p.default_cm,
                save=saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/Dye', saveFile = 'dye_')

    if p.ani_ca2d==True and p.ions_dict['Ca'] == 1 and animate == 1:

        tCa = [1e6*arr[sim.iCa] for arr in sim.cc_time]

        if p.showCells == True:
            viz.AnimateCellData(sim,cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                clrAutoscale = p.autoscale_Ca_ani, clrMin = p.Ca_ani_min_clr, clrMax = p.Ca_ani_max_clr, clrmap = p.default_cm,
                ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Ca',
                saveFile = 'ca_',ignore_simECM = True)
        else:
            viz.AnimateCellData_smoothed(sim,cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                clrAutoscale = p.autoscale_Ca_ani, clrMin = p.Ca_ani_min_clr, clrMax = p.Ca_ani_max_clr, clrmap = p.default_cm,
                ani_repeat=True,number_cells=False,saveFolder = '/animation/Ca', saveFile = 'ca_')

    if p.ani_vm2d==True and animate == 1:

        vmplt = [1000*arr for arr in sim.vm_time]

        if p.sim_ECM == True:

            viz.AnimateCellData(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr,
                clrmap = p.default_cm, ani_repeat=True,number_cells=p.enumerate_cells, current_overlay=p.I_overlay,
                saveFolder = '/animation/Vmem', saveFile = 'vm_')

        elif p.sim_ECM == False:

            if p.showCells == True:
                viz.AnimateCellData(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                     clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr, clrmap = p.default_cm,
                    ani_repeat=True,number_cells=p.enumerate_cells, current_overlay=p.I_overlay,
                    saveFolder = '/animation/Vmem', saveFile = 'vm_')
            else:
                viz.AnimateCellData_smoothed(sim,cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                     clrAutoscale = p.autoscale_Vmem_ani, clrMin = p.Vmem_ani_min_clr, clrMax = p.Vmem_ani_max_clr, clrmap = p.default_cm,
                    ani_repeat=True,number_cells=False,saveFolder = '/animation/Vmem', saveFile = 'vm_',current_overlay=p.I_overlay)

    if p.ani_vmgj2d == True and animate == 1:

        if p.sim_ECM == True:
            viz.AnimateGJData(cells, sim, p, tit='Vcell ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj_', number_cells=False)

        elif p.sim_ECM == False:

            if p.showCells == True:
                viz.AnimateGJData(cells, sim, p, tit='Vmem ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj_', number_cells=False)
            else:
                viz.AnimateGJData_smoothed(cells, sim, p, tit='Vmem ', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                    clrAutoscale = p.autoscale_Vgj_ani, clrMin = p.Vgj_ani_min_clr, clrMax = p.Vgj_ani_max_clr, clrmap = p.default_cm,
                    saveFile = 'vmem_gj', number_cells=False)

    if p.ani_vcell == True and animate == 1 and p.sim_ECM == 1:

        vcellplt = [1000*arr for arr in sim.vcell_time]

        if p.showCells == True:

            viz.AnimateCellData(sim,cells,vcellplt,sim.time,p,tit='V in cell', cbtit = 'Voltage [mV]',
                clrAutoscale = p.autoscale_vcell_ani, clrMin = p.vcell_ani_min_clr, clrMax = p.vcell_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/vcell',
                saveFile = 'vcell_', ignore_simECM =True, current_overlay=p.I_overlay)
        else:
            viz.AnimateCellData_smoothed(sim,cells,vcellplt,sim.time,p,tit='V in cell', cbtit = 'Voltage [mV]',
                clrAutoscale = p.autoscale_vcell_ani, clrMin = p.vcell_ani_min_clr, clrMax = p.vcell_ani_max_clr, clrmap = p.default_cm,
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/vcell', saveFile = 'vcell_',
                current_overlay=p.I_overlay)

    if p.ani_I == True and animate == 1:

        viz.AnimateCurrent(sim,cells,time,p,save=saveAni,ani_repeat=True,current_overlay=p.I_overlay, gj_current =True,
            clrAutoscale=p.autoscale_I_ani,clrMin = p.I_ani_min_clr,clrMax = p.I_ani_max_clr,
            clrmap = p.default_cm, number_cells=False,saveFolder = '/animation/current_gj',saveFile = 'I_')

        if p.sim_ECM == True:

            viz.AnimateCurrent(sim,cells,time,p,save=saveAni,ani_repeat=True,current_overlay=p.I_overlay, gj_current =False,
            clrAutoscale=p.autoscale_I_ani,clrMin = p.I_ani_min_clr,clrMax = p.I_ani_max_clr,
            clrmap = p.default_cm, number_cells=False,saveFolder = '/animation/current_ecm',saveFile = 'I_')


    if p.exportData == True:
        viz.exportData(cells, sim, p)





# def plots4Init(plot_cell,cells,sim,p,saveImages=False):
#
#     if p.plot_single_cell_graphs == True:
#
#         figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,plot_cell,fig=None,
#              ax=None,lncolor='g',ionname='Na+')
#
#         figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,plot_cell,fig=figConcs,
#             ax=axConcs,lncolor='b',ionname='K+')
#
#         figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,plot_cell,fig=figConcs,
#              ax=axConcs,lncolor='r',ionname='M-')
#
#         lg = axConcs.legend()
#         lg.draw_frame(True)
#         titC = 'Concentration of main ions in cell index ' + str(plot_cell) + ' cytoplasm as a function of time'
#         axConcs.set_title(titC)
#         plt.show(block=False)
#
#         figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,plot_cell,fig=None,ax=None,lncolor='b')
#         titV = 'Voltage (Vmem) in cell index ' + str(plot_cell) + ' as a function of time'
#         axVt.set_title(titV)
#         plt.show(block=False)
#
#         if p.ions_dict['Ca'] ==1:
#             figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,plot_cell,fig=None,
#                  ax=None,lncolor='g',ionname='Ca2+ cell')
#             titCa =  'Calcium concentration in cell index ' + str(plot_cell) + ' cytoplasm as a function of time'
#             axA.set_title(titCa)
#             plt.show(block=False)
#
#             if p.Ca_dyn == 1:
#                 figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,plot_cell,fig=None,
#                      ax=None,lncolor='b',ionname='Ca2+ cell')
#                 titER =  'Calcium concentration in cell index ' + str(plot_cell) + ' ER as a function of time'
#                 axD.set_title(titER)
#                 plt.show(block=False)
#
#     if p.plot_vm2d == True:
#
#         if p.sim_ECM == True:
#
#             figV, axV, cbV = viz.plotHetMem(cells,p,zdata=1000*sim.vm_time[-1],number_cells=p.enumerate_cells,
#                 clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm,
#                 edgeOverlay = p.showCells, number_ecm = p.enumerate_cells)
#
#         elif p.sim_ECM == False:
#
#             if p.showCells == True:
#                 figV, axV, cbV = viz.plotPolyData(cells,p,zdata=1000*sim.vm_time[-1],number_cells=p.enumerate_cells,
#                 clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm)
#             else:
#                 figV, axV, cbV = viz.plotCellData(cells,p,zdata=1000*sim.vm_time[-1], clrAutoscale = p.autoscale_Vmem,
#                     clrMin = p.Vmem_min_clr, clrMax = p.Vmem_max_clr, clrmap = p.default_cm, number_cells=p.enumerate_cells)
#
#         axV.set_title('Final Vmem in cell collection')
#         axV.set_xlabel('Spatial distance [um]')
#         axV.set_ylabel('Spatial distance [um]')
#         cbV.set_label('Voltage mV')