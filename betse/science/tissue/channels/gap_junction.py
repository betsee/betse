#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np


class Gap_Junction(object):
    """
    Defines functions controlling gap junction voltage gating characteristics.

    From Harris et al. J of Neurosci.(1983) 3:79-100.
    For Ambystoma Mexicanum (Axolotl!) early embryo gap junction gating.

    """

    def __init__(self, sim, cells, p):

        self.init(sim, cells, p)


    def init(self, sim, cells, p):

        """
        Initialize the gap junction object to its initial conditions

        """

        # parameters from Harris et al 1983
        self.lamb = 0.0013
        self.A1 = 0.077
        self.A2 = 0.14
        self.Ao = 0.217

        # get updated values of the helper-functions:
        # trans-junctional voltage cell 1 wrt 2 to mV units
        V1 = np.zeros(sim.mdl)

        # trans-junctional voltage cell 2 wrt 1 to mV units
        V2 = np.zeros(sim.mdl)

        # voltage-dependent time-constant for channel transitioning from closed to open (1/ms, V in mV):
        self.alpha1 = self.lamb*np.exp(-self.A1*(V1 - p.gj_vthresh))
        self.alpha2 = self.lamb*np.exp(-self.A1*(V2 - p.gj_vthresh))

        # voltage-dependent time constant for channel transitioning from open to closed (1/ms, V in mV):
        beta1 = self.lamb*np.exp(self.A2*(V1 - p.gj_vthresh))
        beta2 = self.lamb*np.exp(self.A2*(V2 - p.gj_vthresh))

        self.beta1 = beta1/(1+50*beta1)
        self.beta2 = beta2/(1+50*beta2)

        # set the initial state to the steady-state conductivity ginf:
        self.n1_inf = self.alpha1/(self.alpha1 + self.beta1)
        self.n2_inf = self.alpha2/ (self.alpha2 + self.beta2)

        self.n1 = self.n1_inf
        self.n2 = self.n2_inf

        gj1 =  self.n1 + p.gj_min*(1 - self.n1)
        gj2 =  self.n2 + p.gj_min*(1 - self.n2)

        gjr = (1/gj1) + (1/gj2)

        sim.gjopen = sim.gj_block*(1/gjr)

    def run(self, sim, cells, p):

        # get updated values of the helper-functions:
        # trans-junctional voltage cell 1 wrt 2 to mV units
        # V1 = (sim.v_cell[cells.nn_i] - sim.v_cell[cells.mem_i])*1e3
        #
        # # # trans-junctional voltage cell 2 wrt 1 to mV units
        # V2 = (sim.v_cell[cells.mem_i] - sim.v_cell[cells.nn_i])*1e3

        V1 = (sim.vm[cells.nn_i] - sim.vm[cells.mem_i])*1e3

        # # trans-junctional voltage cell 2 wrt 1 to mV units
        V2 = (sim.vm[cells.mem_i] - sim.vm[cells.nn_i])*1e3

        # voltage-dependent time-constant for channel transitioning from closed to open (1/ms, V in mV):
        self.alpha1 = self.lamb*np.exp(-self.A1*(V1 - p.gj_vthresh))
        self.alpha2 = self.lamb*np.exp(-self.A1*(V2 - p.gj_vthresh))

        # voltage-dependent time constant for channel transitioning from open to closed (1/ms, V in mV):
        beta1 = self.lamb*np.exp(self.A2*(V1 - p.gj_vthresh))
        beta2 = self.lamb*np.exp(self.A2*(V2 - p.gj_vthresh))

        self.beta1 = beta1/(1+50*beta1)
        self.beta2 = beta2/(1+50*beta2)

        self.delta_n1 = self.alpha1*(1 - self.n1) - (self.n1 - p.gj_min)*self.beta1
        self.delta_n2 = self.alpha2*(1 - self.n2) - (self.n2 - p.gj_min)*self.beta2

        self.n1 = self.n1 + self.delta_n1*p.dt*1e3
        self.n2 = self.n2 + self.delta_n2*p.dt*1e3

        # threshold to ensure 0 to 1 status
        inds_gj_over1 = (self.n1 > 1.0).nonzero()
        self.n1[inds_gj_over1] = 1.0

        inds_gj_under1 = (self.n1 < 0.0).nonzero()
        self.n1[inds_gj_under1] = 0.0

        inds_gj_over2 = (self.n2 > 1.0).nonzero()
        self.n1[inds_gj_over2] = 1.0

        inds_gj_under2 = (self.n2 < 0.0).nonzero()
        self.n2[inds_gj_under2] = 0.0

        # calculate the total channel conductance assuming two independent gj channels, and a gj leak (p.gj_min):
        gj1 = self.n1 + p.gj_min * (1 - self.n1)
        gj2 = self.n2 + p.gj_min * (1 - self.n2)

        # calculate the total (reciprocal) conductance of the channel in terms of the two reciprocal conductances:
        gjr = (1 / gj1) + (1 / gj2)

        # final conductance of the gap junction:
        sim.gjopen = sim.gj_block * (1 / gjr)

        #
        # sim.gjopen = sim.gj_block*(self.n2*self.n1 + p.gj_min*(1 - self.n1)*(1 - self.n2))















