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
import matplotlib.cm as cm
from betse.science import filehandling as fh
import numpy as np
from matplotlib import animation

start_time = time.time()  # get a start value for timing the simulation

#cells = World(vorclose='circle',worldtype='full')  # always need instance of world
#cells.makeWorld()     # call functions to create the world

#fh.saveSim(self.cells.savedWorld,self.cells)   # save the world to cache

#cells = fh.loadWorld(cells.savedWorld)   # load a previously defined world from cache

p = Parameters()

sim = Simulator(p)   # whether running from scratch or loading, instance needs to be called

#sim.baseInit(cells, p)   # initialize data if working from scratch

#sim.runInit(cells,p)     # run and save an initialization if working from scratch

sim,cells, _ = fh.loadSim(sim.savedInit)  # load an initialization from cache

sim.runSim(cells,p,save=False)   # run and save the simulation

#sim,cells,p = fh.loadSim(sim.savedSim)  # load the simulation from cache

vdata_t = np.multiply(sim.vm_time,1000)
# PLOTTING SINGLE CELL DATA
# figC, axC = viz.plotSingleCellCData(self.sim.envcc_time,self.sim.time,self.sim.iNa,self.p.target_cell,fig=None,
#      ax=None,lncolor='g',ionname='Na+')
# figC, axC = viz.plotSingleCellCData(self.sim.envcc_time,self.sim.time,self.sim.iK,self.p.target_cell,fig=figC,
#     ax=axC,lncolor='b',ionname='K+')
# figC, axC = viz.plotSingleCellCData(self.sim.envcc_time,self.sim.time,self.sim.iM,self.p.target_cell,fig=figC,
#      ax=axC,lncolor='r',ionname='M-')

figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iNa,0,fig=None,
     ax=None,lncolor='g',ionname='Na+')
figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iK,0,fig=figC,
    ax=axC,lncolor='b',ionname='K+')
figC, axC = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iM,0,fig=figC,
     ax=axC,lncolor='r',ionname='M-')
lg = axC.legend()
lg.draw_frame(True)
plt.show(block=False)

figCa, axCa = viz.plotSingleCellCData(sim.cc_time,sim.time,sim.iCa,0,
    lncolor='r',ionname='Ca2+')

# plt.figure()
# plt.plot(self.sim.time,self.sim.active_Na_time)
# plt.title('active Na')
# plt.show(block=False)

figVt, axVt = viz.plotSingleCellVData(sim.vm_time,sim.time,0,fig=None,ax=None,lncolor='b')
plt.show(block=False)

# ANIMATING DATA

#viz.AnimateCellData(cells,vdata_t,sim.time,p, save=False, ani_repeat=True,colormap=cm.Blues)

viz.AnimateGJData(cells, sim, p, save=False, ani_repeat=True)
