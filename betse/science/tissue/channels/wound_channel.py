#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Wound-induced transient channel response classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np

from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb


# .................... BASE                               ....................
class WoundABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all wound-induced (TRP based) channel classes.

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
        Initialize wound-induced (TRP based) channel at the point of wounding.

        '''

        self.modulator = 1.0

        V = sim.vm[dyna.targets_vgWound] * 1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)

        self.W_factor = 2.0   # initial concentration of the 'wound factor' (which is basically pressure)
        self.W_decay = 0.1    # FIXME p.W_decay


    def run(self, dyna, sim, cells, p):
        '''
        Simulate wound-induced (TRP based) channel activity after wounding.

        '''

        V = sim.vm[dyna.targets_vgFun] * 1000

        self._calculate_state(V, dyna, sim, p)

        self.W_factor = self.W_factor - self.W_factor*self.W_decay*p.dt

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):
        # calculate m and h channel states using RK4:
        dmWound = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhWound = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Wound = dmWound(dyna.m_Wound, p.dt * 1e3) + dyna.m_Wound
        dyna.h_Wound = dhWound(dyna.h_Wound, p.dt * 1e3) + dyna.h_Wound

        # calculate the open-probability of the channel:
        P = (dyna.m_Wound ** self._mpower) * (dyna.h_Wound ** self._hpower)

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):

        delta_Q = - (dyna.maxDmWound*P*(V - self.vrev))*(self.W_factor/(1 + self.W_factor))

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q, dyna.targets_vgWound, sim, cells, p)
        # self.update_charge(sim.iCa, delta_Q, dyna.targets_vgWound, sim, cells, p)
        # FIXME later on make it optinal to be Na only or  Na + Ca, etc


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
class TRP(WoundABC):
    '''
    TRP model default.

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the wound-induced channel: TRP')

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        dyna.m_Wound = 1.0
        dyna.h_Wound = 1.0

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = 50  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0