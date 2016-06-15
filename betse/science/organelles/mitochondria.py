#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Creates a mitochondria class, which includes a suite of mitochondria-specific molecules, along with
mitochondria-specific pumps (i.e. electron transport chain), channels, and specific methods.
This class also contains the facilities to initialize, define the core computations for a simulation loop,
remove mitochondria during a cutting event, save and report on data, and plot.

"""

import os
import os.path
import numpy as np
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.util.io.log import logs
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionParameters
from betse.science.plot import plot as viz
from betse.science.plot.anim.anim import AnimCellsTimeSeries, AnimEnvTimeSeries
import copy

class Mito(object):

    def __init__(self, sim, cells, p):

        # init basic fields
        self.mit_vol = 0.1*cells.cell_vol     # mit volume
        self.mit_sa = 0.1*cells.cell_sa      # mit surface areas
        self.Vmit = np.zeros(sim.cdl)   # transmembrane voltage for mit
        self.Q = np.zeros(sim.cdl)     # total charge in mit
        self.cm_mit = self.mit_sa*p.cm    # mit membrane capacitance

        sim.cc_mit = copy.deepcopy(sim.cc_cells)    # ion concentrations
        self.Dm_mit = copy.deepcopy(sim.cc_cells)    # membrane permeability

        for arr in self.Dm_mit:

            arr[:] = 1.0e-18                 # membrane permeability altered so all are minimal

        self.Dm_mit_base = copy.deepcopy(self.Dm_mit)  # copies of Dm for ion channel dynamics
        self.Dm_channels = copy.deepcopy(self.Dm_mit)


    def get_v(self, sim, p):

        self.Q = stb.get_charge(sim.cc_mit, sim.zs, self.mit_vol, p)
        self.Vmit = (1/self.cm_mit)*self.Q

    def update(self, sim, cells, p): # FIXME electroflux between cell and env for all ions, get_v


        # self.channels(sim, cells, p)

        # flux = stb.electroflux()
        # update with flux

        # self.get_v(sim, p)

        pass

    def channels(self, sim, cells, p):

        pass   # FIXME finish this







