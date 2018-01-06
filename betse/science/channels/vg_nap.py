#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated sodium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channels.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

# .................... BASE                               ....................
class VgNaPABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all persistent Voltage-gated sodium channel classes.

    Attributes
    ----------
    _mInf :
        Equillibrium value function for m-gates of channel
    _mTau :
        Time-constant function for m-gates of channel
    _hInf :
        Equillibrium value function for h-gates of channel
    _hTau :
        Time-constant function for h-gates of channel
    '''

    def init(self, dyna, sim, cells, p):
        '''
        Initialize targeted voltage-gated sodium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        self.modulator = 1.0

        self.v_corr = 0.0  # offset of voltages in the model -- experimental junction voltage [mV]

        V = sim.vm*1000 + self.v_corr

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''
        V = sim.vm* 1000 + self.v_corr

        self._calculate_state(V, dyna, sim, p)

        self._implement_state(V, dyna,sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):

        # calculate m and h channel states using RK4:
        dmNaP = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhNaP = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_NaP = dmNaP(dyna.m_NaP, p.dt*1e3) + dyna.m_NaP
        dyna.h_NaP = dhNaP(dyna.h_NaP, p.dt*1e3) + dyna.h_NaP

        # calculate the open-probability of the channel:
        P = (dyna.m_NaP ** self._mpower) * (dyna.h_NaP ** self._hpower)

        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the original Hodgkin Huxley equation.

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        #delta_Q = -(dyna.maxDmNaP * P * (V - self.vrev))

        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem = sim.cc_cells[sim.iNa][cells.mem_to_cells]

        if p.is_ecm is True:
            c_env = sim.cc_env[sim.iNa][cells.map_mem2ecm]

        else:
            c_env = sim.cc_env[sim.iNa]

        IdM = np.ones(sim.mdl)

        z_ion = sim.zs[sim.iNa] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmNaP*P*self.modulator*sim.rho_channel

        # calculate specific ion flux contribution for this channel:
        delta_Q = stb.electroflux(c_env, c_mem, Dchan, p.tm * IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.chan_flux = np.zeros(sim.mdl)
        self.chan_flux[dyna.targets_vgNaP] = -delta_Q[dyna.targets_vgNaP]

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q, dyna.targets_vgNaP, sim, cells, p)


    @abstractmethod
    def _init_state(self, V, dyna, sim, p):
        '''
        Do something.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V, dyna, sim, p):
        '''
        Do something.
        '''
        pass

# ....................{ SUBCLASS                           }....................

class Nav1p6(VgNaPABC):

    """
    Nav1.6 was detected during the embryonic period in brain. This is the most abundantly expressed
    NaV channel in the CNS during adulthood. Is the most abundant channel at mature nodes of Ranvier
    in myelinated axons in the CNS. NaV1.6 is broadly expressed in the nervous system in a variety
    of cells including Purkinje cells, motor neurons, pyramidal and granule neurons, glial cells and
    Schwann cells and is enriched at the nodes of Ranvier.
    Nav1.6 channels have been also detected in immune cells, such microglia and macrophagues
     and in cultured microglia, Nav1.6 is the most prominently expressed sodium channel.

     Reference: Smith MR. et al. Functional analysis of the mouse Scn8a sodium channel. \
     J. Neurosci., 1998 Aug 15 , 18 (6093-102).

    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the persistent vgNa channel: Nav1p6')

        self.vrev = 50     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        self._mpower = 1.0
        self._hpower = 0.0

        dyna.m_NaP = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        dyna.h_NaP = 1.0

    def _calculate_state(self, V, dyna, sim, p):

        self._mInf = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self._mTau = 1

        self._hInf = 1.0
        self._hTau = 1.0