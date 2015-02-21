#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import matplotlib.pyplot as plt
import time
from betse.science.world import World
from betse.science.compute import Simulator
from betse.science.parameters import Parameters
from betse.science import visualize as viz
from betse.science import filehandling as fh
import matplotlib.cm as cm
import numpy as np
import os, os.path


class SimRunner(object):

    def initialize(self):

        #.....COMMAND SEQUENCE #1: "initialize" ...........................................................................
        # Run an initialization simulation from scratch and save it to the initialization cache.

        start_time = time.time()  # get a start value for timing the simulation

        cells = World(vorclose='circle',worldtype='full')  # create an instance of world
        cells.makeWorld()     # call function to create the world
        p = Parameters()     # create an instance of Parameters
        p.time_profile = 'initialize'  # enforce the time profile to be initialize
        sim = Simulator(p)   # create an instance of Simulator
        sim.baseInit(cells, p)   # initialize simulation data structures
        sim.runInit(cells,p)     # run and save the initialization

        print('The initialization took', round(time.time() - start_time,2), 'seconds to complete')

        plots4Init(p.plot_cell,cells,sim,p,saveImages=p.autosave)

        plt.show()

    def simulate(self):

        # COMMAND SEQUENCE #2: "runSim"  ...............................................................................
        # Run simulation from a previously saved initialization. FIXME: throw exception if init cache is empty.
        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters()     # create an instance of Parameters
        p.time_profile = 'simulate'  # enforce the time-profile to be simulate
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        sim.runSim(cells,p,save=True)   # run and optionally save the simulation to the cache

        print('The simulation took', round(time.time() - start_time,2), 'seconds to complete')

        plots4Sim(p.plot_cell,cells,sim,p,saveImages = p.autosave,animate=p.createAnimations,saveAni=p.saveAnimations)

        plt.show()

    def loadInit(self):
        # COMMAND SEQUENCE #3: "loadInit" ..............................................................................
        # Load and visualize a previously solved initialization

        p = Parameters()     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

        plots4Init(p.plot_cell,cells,sim,p,saveImages=p.autosave)

        plt.show()

    def loadSim(self):

        # COMMAND SEQUENCE #4: "loadSim" ...............................................................................
        # Load and visualize a previously solved simulation

        p = Parameters()     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells,_ = fh.loadSim(sim.savedSim)  # load the simulation from cache

        plots4Sim(p.plot_cell,cells,sim,p,saveImages=p.autosave, animate=p.createAnimations, saveAni=p.saveAnimations)

        plt.show()

def plots4Init(plot_cell,cells,sim,p,saveImages=False):

    if p.plot_single_cell_graphs == True:

        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,plot_cell,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,plot_cell,fig=figConcs,
            ax=axConcs,lncolor='b',ionname='K+')
        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,plot_cell,fig=figConcs,
             ax=axConcs,lncolor='r',ionname='M-')
        lg = axConcs.legend()
        lg.draw_frame(True)
        titC = 'Concentration of main ions in cell index ' + str(plot_cell) + ' cytoplasm as a function of time'
        axConcs.set_title(titC)
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,plot_cell,fig=None,ax=None,lncolor='b')
        titV = 'Voltage (Vmem) in cell index ' + str(plot_cell) + ' as a function of time'
        axVt.set_title(titV)
        plt.show(block=False)

        if p.ions_dict['Ca'] ==1:

            figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,plot_cell,fig=None,
                 ax=None,lncolor='g',ionname='Ca2+ cell')
            titCa =  'Calcium concentration in cell index ' + str(plot_cell) + ' cytoplasm as a function of time'
            axA.set_title(titCa)
            plt.show(block=False)

            if p.Ca_dyn == 1:

                figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,plot_cell,fig=None,
                     ax=None,lncolor='b',ionname='Ca2+ cell')
                titER =  'Calcium concentration in cell index ' + str(plot_cell) + ' ER as a function of time'
                axD.set_title(titER)
                plt.show(block=False)

    if p.plot_vm2d == True:
        if p.showCells == True:
            figV, axV, cbV = viz.plotPolyData(cells,p,zdata=1000*sim.vm_time[-1],number_cells=p.enumerate_cells)
        else:
            figV, axV, cbV = viz.plotCellData(cells,p,zdata=1000*sim.vm_time[-1])

        axV.set_title('Final Vmem in cell collection')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')


def plots4Sim(plot_cell,cells,sim,p, saveImages=False, animate=0,saveAni=False):

    if saveImages == True:
        images_path = p.sim_results
        image_cache_dir = os.path.expanduser(images_path)
        os.makedirs(image_cache_dir, exist_ok=True)
        savedImg = os.path.join(image_cache_dir, 'fig_')

    if p.plot_single_cell_graphs == True:

        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,plot_cell,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,plot_cell,fig=figConcs,
            ax=axConcs,lncolor='b',ionname='K+')
        figConcs, axConcs = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,plot_cell,fig=figConcs,
             ax=axConcs,lncolor='r',ionname='M-')
        lg = axConcs.legend()
        lg.draw_frame(True)
        titC = 'Main ions in cell index ' + str(plot_cell)
        axConcs.set_title(titC)
        if saveImages == True:
            savename1 = savedImg + 'conc_time'
            plt.savefig(savename1,dpi=300,format='png')
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,plot_cell,fig=None,ax=None,lncolor='b')
        titV = 'Voltage (Vmem) in cell ' + str(plot_cell)
        axVt.set_title(titV)
        if saveImages == True:
            savename2 = savedImg + 'Vmem_time'
            plt.savefig(savename2,dpi=300,format='png')
        plt.show(block=False)

        if p.ions_dict['Ca'] ==1:

            figA, axA = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,plot_cell,fig=None,
                 ax=None,lncolor='g',ionname='Ca2+ cell')
            titCa =  'Cytosolic Ca2+ in cell index ' + str(plot_cell)
            axA.set_title(titCa)

            if saveImages == True:
                savename3 = savedImg + 'cytosol_Ca_time'
                plt.savefig(savename3,dpi=300,format='png')

            plt.show(block=False)

            if p.Ca_dyn == 1:

                figD, axD = viz.plotSingleCellCData(sim.cc_er_time,sim.time,0,plot_cell,fig=None,
                     ax=None,lncolor='b',ionname='Ca2+ cell')
                titER =  'ER Ca2+ in cell index ' + str(plot_cell)
                axD.set_title(titER)

                if saveImages == True:
                    savename4 = savedImg + 'ER_Ca_time'
                    plt.savefig(savename4,dpi=300,format='png')

                plt.show(block=False)

    if p.plot_vm2d == True:

        if p.showCells == True:
            figV, axV, cbV = viz.plotPolyData(cells,p,zdata=1000*sim.vm_time[-1],number_cells=p.enumerate_cells)
        else:
            figV, axV, cbV = viz.plotCellData(cells,p,zdata=1000*sim.vm_time[-1])

        axV.set_title('Final Vmem')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        if saveImages == True:
            savename5 = savedImg + 'final_Vmem_2D'
            plt.savefig(savename5,dpi=300,format='png')

        plt.show(block=False)


    if  p.plot_ip32d == True and p.scheduled_options['IP3'] != 0:

        if p.showCells == True:
            figIP3, axIP3, cbIP3 = viz.plotPolyData(cells,p,zdata=sim.cIP3_time[-1]*1e3,number_cells=p.enumerate_cells)
        else:
             figIP3, axIP3, cbIP3 = viz.plotCellData(cells,p,zdata=sim.cIP3_time[-1]*1e3)

        axIP3.set_title('Final IP3 concentration')
        axIP3.set_xlabel('Spatial distance [um]')
        axIP3.set_ylabel('Spatial distance [um]')
        cbIP3.set_label('Concentration umol/L')

        if saveImages == True:
            savename6 = savedImg + 'final_IP3_2D'
            plt.savefig(savename6,dpi=300,format='png')

        plt.show(block=False)

    if p.plot_dye2d == True and p.voltage_dye == 1:

        if p.showCells == True:
            figVdye, axVdye, cbVdye = viz.plotPolyData(cells,p,zdata=sim.cDye_time[-1]*1e3,number_cells=p.enumerate_cells)
        else:
            figVdye, axVdye, cbVdye = viz.plotCellData(cells,p,zdata=sim.cDye_time[-1]*1e3)

        axVdye.set_title('Final voltage-sensitive dye')
        axVdye.set_xlabel('Spatial distance [um]')
        axVdye.set_ylabel('Spatial distance [um]')
        cbVdye.set_label('Concentration umol/L')

        if saveImages == True:
            savename7 = savedImg + 'final_dye_2D'
            plt.savefig(savename7,dpi=300,format='png')

        plt.show(block=False)

    if p.plot_ca2d ==True and  p.ions_dict['Ca'] == 1:

        if p.showCells == True:
            figCa, axCa, cbCa = viz.plotPolyData(cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6,number_cells= p.enumerate_cells)
        else:
            figCa, axCa, cbCa = viz.plotCellData(cells,p,zdata=sim.cc_time[-1][sim.iCa]*1e6)

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if saveImages == True:
            savename8 = savedImg + 'final_Ca_2D'
            plt.savefig(savename8,dpi=300,format='png')

        plt.show(block=False)


    if p.ani_ip32d ==True and p.scheduled_options['IP3'] != 0 and animate == 1:

        IP3plotting = np.asarray(sim.cIP3_time)
        IP3plotting = np.multiply(IP3plotting,1e3)

        if p.showCells == True:
            viz.AnimateCellData(cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                save= saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/IP3', saveFile = 'ip3_')

        else:
            viz.AnimateCellData_smoothed(cells,IP3plotting,sim.time,p,tit='IP3 concentration', cbtit = 'Concentration [umol/L]',
                save= saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/IP3', saveFile = 'ip3_')

    if p.ani_dye2d == True and p.voltage_dye == 1 and animate ==1:

        Dyeplotting = np.asarray(sim.cDye_time)
        Dyeplotting = np.multiply(Dyeplotting,1e3)

        if p.showCells == True:
            viz.AnimateCellData(cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                save=saveAni, ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Dye', saveFile = 'dye_')

        else:
            viz.AnimateCellData_smoothed(cells,Dyeplotting,sim.time,p,tit='V-sensitive dye', cbtit = 'Concentration [umol/L]',
                save=saveAni, ani_repeat=True,number_cells=False,saveFolder = '/animation/Dye', saveFile = 'dye_')


    if p.ani_ca2d==True and p.ions_dict['Ca'] == 1 and animate == 1:

        tCa = [1e6*arr[sim.iCa] for arr in sim.cc_time]

        if p.showCells == True:
            viz.AnimateCellData(cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Ca', saveFile = 'ca_')
        else:
            viz.AnimateCellData_smoothed(cells,tCa,sim.time,p,tit='Cytosolic Ca2+', cbtit = 'Concentration [nmol/L]', save=saveAni,
                ani_repeat=True,number_cells=False,saveFolder = '/animation/Ca', saveFile = 'ca_')

    if p.ani_vm2d==True and animate == 1:

        vmplt = [1000*arr for arr in sim.vm_time]

        if p.showCells == True:
            viz.AnimateCellData(cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                ani_repeat=True,number_cells=p.enumerate_cells,saveFolder = '/animation/Vmem', saveFile = 'vm_')
        else:
            viz.AnimateCellData_smoothed(cells,vmplt,sim.time,p,tit='Cell Vmem', cbtit = 'Voltage [mV]', save=saveAni,
                ani_repeat=True,number_cells=False,saveFolder = '/animation/Vmem', saveFile = 'vm_')

    if p.ani_vmgj2d == True and animate == 1:

        if p.showCells == True:
            viz.AnimateGJData(cells, sim, p, tit='Cell Vmem', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                saveFile = 'vmem_gj_', number_cells=False)

        else:
            viz.AnimateGJData_smoothed(cells, sim, p, tit='Cell Vmem', save=saveAni, ani_repeat=True,saveFolder = '/animation/Vmem_gj',
                saveFile = 'vmem_gj', number_cells=False)

    if p.exportData == True:
        viz.exportData(cells, sim, p)

boo = SimRunner()
boo.simulate()

