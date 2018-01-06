#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Creates a mitochondria class, which includes a suite of mitochondria-specific molecules, along with
mitochondria-specific pumps (i.e. electron transport chain), channels, and specific methods.
This class also contains the facilities to initialize, define the core computations for a simulation loop,
remove mitochondria during a cutting event, save and report on data, and plot.

"""

import copy

import numpy as np

from betse.science import sim_toolbox as stb


class Mito(object):

    def __init__(self, sim, cells, p):

        # init basic fields
        self.mit_vol = 0.5*cells.cell_vol     # mit volume
        # print("!!!!!!!!! mit_vol: {}".format(self.mit_vol))
        self.mit_sa = 0.5*cells.cell_sa      # mit surface areas
        # self.Vmit = np.zeros(sim.cdl)   # initial transmembrane voltage for mit
        self.Vmit = np.zeros(sim.cdl)  # initial transmembrane voltage for mit
        self.Q = np.zeros(sim.cdl)     # total charge in mit
        self.cm_mit = self.mit_sa*p.cm    # mit membrane capacitance

        sim.cc_mit = copy.deepcopy(sim.cc_cells)    # ion concentrations
        self.Dm_mit = copy.deepcopy(sim.cc_cells)    # membrane permeability

        for arr in self.Dm_mit:

            arr[:] = 1.0e-19                 # membrane permeability altered so all are minimal

        # set calcium concentration in mitochondria to an initially low value:
        # if p.ions_dict['Ca'] == 1:
        #     self.Dm_mit[sim.iCa] = 1.0e-15  # add a mitochondrial calcium uniporter set
        #     # sim.cc_mit[sim.iCa] = 10.0e-6  # [100 nM]

        self.Dm_mit_base = copy.deepcopy(self.Dm_mit)  # copies of Dm for ion channel dynamics
        self.Dm_channels = copy.deepcopy(self.Dm_mit)

        self.zmit = copy.deepcopy(sim.cc_cells)

        for i, arr in enumerate(self.zmit):
            arr[:] = sim.zs[i]

    def get_v(self, sim, p):

        self.Q = stb.get_charge(sim.cc_mit, self.zmit, self.mit_vol, p) + self.extra_rho*self.mit_vol

        self.Vmit = (1/self.cm_mit)*self.Q

    def update(self, sim, cells, p):

        # self.channels(sim, cells, p)

        for i in sim.movingIons:

            IdCM = np.ones(sim.cdl)

            f_ED = stb.electroflux(sim.cc_cells[i], sim.cc_mit[i], self.Dm_mit[i], p.tm*IdCM, sim.zs[i]*IdCM,
                self.Vmit, sim.T, p, rho=1)

            # update with flux
            sim.cc_cells[i] = sim.cc_cells[i] - f_ED*(self.mit_sa/cells.cell_vol)*p.dt
            sim.cc_mit[i] = sim.cc_mit[i] + f_ED*(self.mit_sa/self.mit_vol)*p.dt

        self.get_v(sim, p)

    # def channels(self, sim, cells, p):

    # FIXME major useful channel would be the Ca2+ uniporter, which opens with 1 - 5 um cytosolic calcium
    #
    #     pass   # FIXME finish this

    def remove_mits(self, sim, target_inds_cell):

        # remove cells from the mit voltage list:
        vmit2 = np.delete(self.Vmit, target_inds_cell)
        # reassign the new data vector to the object:
        self.Vmit = vmit2

        mitv2 = np.delete(self.mit_vol, target_inds_cell)
        self.mit_vol = mitv2

        mitca2 = np.delete(self.mit_sa, target_inds_cell)
        self.mit_sa = mitca2

        Q2 = np.delete(self.Q, target_inds_cell)
        self.Q = Q2

        cm2 = np.delete(self.cm_mit, target_inds_cell)
        self.cm_mit = cm2

        cc_mit2 = []

        for i, arr in enumerate(sim.cc_mit):

            # remove cells from the mit ion array in sim:
            arr2 = np.delete(arr, target_inds_cell)
            cc_mit2.append(arr2)

        sim.cc_mit = np.asarray(cc_mit2)

        dmit1 = []
        dmit2 = []
        dmit3 = []
        zmit = []

        for i, arr in enumerate(self.Dm_mit):

            arr2 = np.delete(arr, target_inds_cell)
            dmit2.append(arr2)

        self.Dm_mit = np.asarray(dmit2)

        for i, arr in enumerate(self.Dm_mit_base):

            arr2 = np.delete(arr, target_inds_cell)
            dmit1.append(arr2)

        self.Dm_mit_base = np.asarray(dmit1)

        for i, arr in enumerate(self.Dm_channels):

            arr2 = np.delete(arr, target_inds_cell)
            dmit3.append(arr2)

        self.Dm_channels = np.asarray(dmit3)

        for i, arr in enumerate(self.zmit):

            arr2 = np.delete(arr, target_inds_cell)
            zmit.append(arr2)

        self.zmit = np.asarray(zmit)







