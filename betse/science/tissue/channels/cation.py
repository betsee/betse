#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Cation leak channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np
from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb


# .................... BASE                               ....................
class CationABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all non-selective cation channel classes.

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
        Initialize targeted voltage-gated funny current channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        self.modulator = 1.0

        V = sim.vm*1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated funny current channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        V = sim.vm*1000

        self._calculate_state(V, dyna, sim, p)

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):
        # calculate m and h channel states using RK4:
        dmCat = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhCat = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Cat = dmCat(dyna.m_Cat, p.dt * 1e3) + dyna.m_Cat
        dyna.h_Cat = dhCat(dyna.h_Cat, p.dt * 1e3) + dyna.h_Cat

        # calculate the open-probability of the channel:
        P = (dyna.m_Cat ** self._mpower) * (dyna.h_Cat ** self._hpower)

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        # delta_Q = - (dyna.maxDmCat*P*(V - self.vrev))

        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem_Na = sim.cc_cells[sim.iNa][cells.mem_to_cells]
        c_mem_K = sim.cc_cells[sim.iK][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env_Na = sim.cc_env[sim.iNa][cells.map_mem2ecm]
            c_env_K = sim.cc_env[sim.iK][cells.map_mem2ecm]

        else:
            c_env_Na = sim.cc_env[sim.iNa]
            c_env_K = sim.cc_env[sim.iK]

        IdM = np.ones(sim.mdl)

        z_Na = sim.zs[sim.iNa] * IdM
        z_K = sim.zs[sim.iK] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmCat * P * 1.0e-9

        # calculate specific ion flux contribution for this channel:
        delta_Q_Na = stb.electroflux(c_env_Na, c_mem_Na, Dchan, p.tm * IdM, z_Na, sim.vm, sim.T, p, rho=sim.rho_channel)
        delta_Q_K = stb.electroflux(c_env_K, c_mem_K, Dchan, p.tm * IdM, z_K, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.clip_flux(delta_Q_Na, threshold=p.flux_threshold)
        self.clip_flux(delta_Q_K, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q_Na, dyna.targets_vgCat, sim, cells, p)
        self.update_charge(sim.iK, delta_Q_K, dyna.targets_vgCat, sim, cells, p)


    @abstractmethod
    def _init_state(self, V, dyna, sim, p):
        '''
        Initializes values of the m and h gates of the channel.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V, dyna, sim, p):
        '''
        Calculates time-dependent values of the m and h gates of the channel.
        '''
        pass

# ....................{ SUBCLASS                           }....................
class Cat_Leak(CationABC):
    '''

    Membrane leak channel letting Na+, K+ and Ca2+ into the cell

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the cation channel: Cat_Leak')

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        dyna.m_Cat = 1.0000 / (1 + np.exp((V - -99) / 6.2))
        dyna.h_Cat = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0