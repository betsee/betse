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

        print('The initialization took', time.time() - start_time, 'seconds to complete')

        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,0,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,0,fig=figC,
            ax=axC,lncolor='b',ionname='K+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,0,fig=figC,
             ax=axC,lncolor='r',ionname='M-')
        lg = axC.legend()
        lg.draw_frame(True)
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,0,fig=None,ax=None,lncolor='b')
        plt.show(block=False)

        plt.show()

    def simulate(self, savePNG=False):

        # COMMAND SEQUENCE #2: "runSim"  ...............................................................................
        # Run simulation from a previously saved initialization. FIXME: throw exception if init cache is empty.
        start_time = time.time()  # get a start value for timing the simulation

        p = Parameters()     # create an instance of Parameters
        p.time_profile = 'simulate'  # enforce the time-profile to be simulate
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache
        sim.runSim(cells,p,save=False)   # run and optionally save the simulation to the cache

        print('The simulation took', time.time() - start_time, 'seconds to complete')

        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,0,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,0,fig=figC,
            ax=axC,lncolor='b',ionname='K+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,0,fig=figC,
             ax=axC,lncolor='r',ionname='M-')
        lg = axC.legend()
        lg.draw_frame(True)
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,0,fig=None,ax=None,lncolor='b')
        plt.show(block=False)

        viz.AnimateGJData(cells, sim, p, save=False, ani_repeat=True)
        #viz.AnimateCellData(cells,sim.active_K_time,sim.time,p, save=False,ani_repeat=True)
        if savePNG == True:
            viz.Animate2PNG(cells,sim,p)

        plt.show()

    def loadInit(self):
        # COMMAND SEQUENCE #3: "loadInit" ..............................................................................
        # Load and visualize a previously solved initialization

        p = Parameters()     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells, _ = fh.loadSim(sim.savedInit)  # load the initialization from cache

        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,0,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,0,fig=figC,
            ax=axC,lncolor='b',ionname='K+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,0,fig=figC,
             ax=axC,lncolor='r',ionname='M-')
        lg = axC.legend()
        lg.draw_frame(True)
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,0,fig=None,ax=None,lncolor='b')
        plt.show(block=False)

        plt.show()

    def loadSim(self):

        # COMMAND SEQUENCE #4: "loadSim" ...............................................................................
        # Load and visualize a previously solved simulation

        p = Parameters()     # create an instance of Parameters
        sim = Simulator(p)   # create an instance of Simulator
        sim,cells,p = fh.loadSim(sim.savedSim)  # load the simulation from cache

        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,0,fig=None,
             ax=None,lncolor='g',ionname='Na+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,0,fig=figC,
            ax=axC,lncolor='b',ionname='K+')
        figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,0,fig=figC,
             ax=axC,lncolor='r',ionname='M-')
        lg = axC.legend()
        lg.draw_frame(True)
        plt.show(block=False)

        figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,0,fig=None,ax=None,lncolor='b')
        plt.show(block=False)

        plt.show()

boo = SimRunner()
boo.simulate()

