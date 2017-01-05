#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated chloride channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np
from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb


class VgClABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated chloride channel classes.

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
        Initialize targeted voltage-gated potassium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        self.modulator = 1.0

        self.v_corr = 0.0   # in experiments, the measurement junction voltage is about 10 mV

        V = sim.vm*1000 + self.v_corr

        self._init_state(V, dyna, sim, p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        V = sim.vm*1000 + self.v_corr

        self._calculate_state(V, dyna, sim, p)

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):

        # calculate m and h channel states using RK4:
        dmCl = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhCl = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Cl = dmCl(dyna.m_Cl, p.dt * 1e3) + dyna.m_Cl
        dyna.h_Cl = dhCl(dyna.h_Cl, p.dt * 1e3) + dyna.h_Cl

        # calculate the open-probability of the channel:
        P = (dyna.m_Cl ** self._mpower) * (dyna.h_Cl ** self._hpower)

        # print(P.min(), P.max(), P.mean())

        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the original Hodgkin Huxley equation.

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        # delta_Q = - (dyna.maxDmK*P*(V - self.vrev))

        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem = sim.cc_cells[sim.iCl][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env = sim.cc_env[sim.iCl][cells.map_mem2ecm]

        else:
            c_env = sim.cc_env[sim.iCl]

        IdM = np.ones(sim.mdl)

        z_ion = sim.zs[sim.iCl] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmCl*P*1.0e-9    # 1.0e-9 multiplier to approximately convert from conductivity

        # calculate specific ion flux contribution for this channel:
        delta_Q = stb.electroflux(c_env, c_mem, Dchan, p.tm * IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iCl, delta_Q, dyna.targets_vgCl, sim, cells, p)


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


class ClLeak(VgClABC):

    '''
    Simple chloride leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a substance-modulated Cl- channel')


        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Cl = np.ones(sim.mdl)
        dyna.h_Cl = np.ones(sim.mdl)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1
        self._mTau = 1
        self._hInf = 1
        self._hTau = 1

