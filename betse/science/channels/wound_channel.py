#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Wound-induced transient channel response classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.exceptions import BetseSimConfException
from betse.science import sim_toolbox as stb
from betse.science.channels.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs
# from betse.science.chemistry.molecule import get_influencers

# .................... BASE                               ....................
class WoundABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all wound-induced (TRP-based) channel classes.

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

    def init(self, vm, cells, p, targets = None):
        '''
        Initialize wound-induced (TRP based) channel at the point of wounding.

        '''

        if targets is None:

            self.targets = cells.mem_i

        else:
            self.targets = targets

        self.data_length = len(self.targets)
        self.mdl = len(cells.mem_i)

        self.modulator = 1.0

        self.v_corr = 0.0   # in experiments, the measurement junction voltage is about 10 mV

        V = vm[self.targets]*1000 + self.v_corr
        # V = vm * 1000 + self.v_corr

        self._init_state(V)

        self.ions = ['Na', 'K']
        self.rel_perm = [1.0, 1.0]

        self.W_factor = 2.0   # initial concentration of the 'wound factor' (which is basically pressure)
        self.W_decay = p.wound_close_factor


    def run(self, vm, p):
        '''
        Simulate wound-induced (TRP based) channel activity after wounding.

        '''

        V = vm[self.targets]*1000

        self._calculate_state(V)

        self.W_factor = self.W_factor - self.W_factor*self.W_decay*p.dt

        self._implement_state(V, p)

    def _implement_state(self, V, p):
        # Update the 'm' and 'h' gating functions:
        self.update_mh(p, time_unit = self.time_unit)

        # calculate the open-probability of the channel:
        P = (self.m ** self._mpower) * (self.h ** self._hpower)

        self.P = np.zeros(self.mdl)
        self.P[self.targets] = P


        # # obtain concentration of ion inside and out of the cell, as well as its charge z:
        # c_mem_Na = sim.cc_cells[sim.iNa][cells.mem_to_cells]
        # c_mem_K = sim.cc_cells[sim.iK][cells.mem_to_cells]
        #
        # if p.is_ecm is True:
        #     c_env_Na = sim.cc_env[sim.iNa][cells.map_mem2ecm]
        #     c_env_K = sim.cc_env[sim.iK][cells.map_mem2ecm]
        #
        # else:
        #     c_env_Na = sim.cc_env[sim.iNa]
        #     c_env_K = sim.cc_env[sim.iK]
        #
        # IdM = np.ones(sim.mdl)
        #
        # z_Na = sim.zs[sim.iNa] * IdM
        # z_K = sim.zs[sim.iK] * IdM
        #
        # # membrane diffusion constant of the channel:
        # Dchan = dyna.maxDmWound*P*(self.W_factor/(1 + self.W_factor))*self.modulator
        #
        # # calculate specific ion flux contribution for this channel:
        # delta_Q_Na = stb.electroflux(c_env_Na, c_mem_Na, Dchan, p.tm * IdM, z_Na, sim.vm, sim.T, p, rho=sim.rho_channel)
        # delta_Q_K = stb.electroflux(c_env_K, c_mem_K, Dchan, p.tm * IdM, z_K, sim.vm, sim.T, p, rho=sim.rho_channel)
        #
        # self.clip_flux(delta_Q_Na, threshold=p.flux_threshold)
        # self.clip_flux(delta_Q_K, threshold=p.flux_threshold)
        #
        # self.chan_flux = np.zeros(sim.mdl)
        # self.chan_flux[dyna.targets_vgWound] = -delta_Q_Na[dyna.targets_vgWound] -delta_Q_K[dyna.targets_vgWound]
        #
        # self.update_charge(sim.iNa, delta_Q_Na, dyna.targets_vgWound, sim, cells, p)
        # self.update_charge(sim.iK, delta_Q_K, dyna.targets_vgWound, sim, cells, p)

        # if p.ions_dict['Ca'] == 1.0:
        #     self.update_charge(sim.iCa, 0.1*delta_Q, dyna.targets_vgWound, sim, cells, p)


    @abstractmethod
    def _init_state(self, V, sim, p):
        '''
        Initializes values of the m and h gates of the channel.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V, sim, p):
        '''
        Calculates time-dependent values of the m and h gates of the channel.
        '''
        pass

# ....................{ SUBCLASS                           }....................
class TRP(WoundABC):
    '''
    TRP model default.

    '''

    def _init_state(self, V, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the wound-induced channel: TRP')

        self.v_corr = 0

        # initialize values of the m and h gates of the channel on m_inf and h_inf:
        self.m = np.ones(sim.mdl)
        self.h = np.ones(sim.mdl)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = 0  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0


