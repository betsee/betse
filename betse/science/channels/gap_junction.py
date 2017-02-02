#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
        # voltage-dependent time-constant for channel transitioning from closed to open (1/ms, V in mV):

        # we use the absolute value of transjunctional voltage due to symmetric channels and convert to mV:
        V1 = 1.0e3*np.abs(sim.vgj)

        # time constant transitioning from closed to open:
        self.alpha = self.lamb*np.exp(-self.A1*(V1 - p.gj_vthresh))

        # voltage-dependent time constant for channel transitioning from open to closed (1/ms, V in mV):
        beta = self.lamb*np.exp(self.A2*(V1 - p.gj_vthresh))

        self.beta = beta/(1+50*beta)

        # set the initial state to the steady-state conductivity ginf:
        self.gjopen = sim.gj_block*((self.alpha/(self.alpha + self.beta)) + p.gj_min)


    def run(self, sim, cells, p):


        # we use the absolute value of transjunctional voltage due to symmetric channels and convert to mV:
        V1 = 1.0e3 * np.abs(sim.vgj)

        # time constant transitioning from closed to open:
        self.alpha = self.lamb * np.exp(-self.A1 * (V1 - p.gj_vthresh))

        # voltage-dependent time constant for channel transitioning from open to closed (1/ms, V in mV):
        beta = self.lamb * np.exp(self.A2 * (V1 - p.gj_vthresh))

        self.beta = beta / (1 + 50 * beta)

        # update GJ conductivity with incremental change of channel state (time is converted to be in ms):
        sim.gjopen = sim.gjopen + ((1.0 - sim.gjopen)*self.alpha - (sim.gjopen - p.gj_min)*self.beta)*p.dt*1e3

        # apply any blocking modulation:
        sim.gjopen = sim.gj_block*sim.gjopen

        # print(sim.gjopen.mean())


class Gap_Junction_i(object):
    """
    Defines simple gap junction voltage gating characteristics.

    """

    def __init__(self, sim, cells, p):

        self.init(sim, cells, p)


    def init(self, sim, cells, p):

        """
        Initialize the gap junction object to its initial conditions

        """

        sim.gjopen = p.gj_min +(-1 / (1 + np.exp((sim.vgj*1e3 + p.gj_vthresh) / p.gj_vgrad)) +
                      1 / (1 + np.exp((sim.vgj*1e3 - p.gj_vthresh) / p.gj_vgrad)))*sim.gj_block

    def run(self, sim, cells, p):


        sim.gjopen = p.gj_min +(-1 / (1 + np.exp((sim.vgj*1e3 + p.gj_vthresh) / p.gj_vgrad)) +
                      1 / (1 + np.exp((sim.vgj*1e3 - p.gj_vthresh) / p.gj_vgrad)))*sim.gj_block















