#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated sodium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np

from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb


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

    def init(self, dyna, sim, p):
        '''
        Initialize targeted voltage-gated sodium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        V = sim.vm[dyna.targets_vgNaP] * 1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        self._calculate_state(
            V=sim.vm[dyna.targets_vgNaP] * 1000,
            dyna=dyna, sim=sim, p=p)

        self._implement_state(dyna,sim,p)

    def _implement_state(self, dyna, sim, p):

        # calculate m and h channel states using RK4:
        dmNaP = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhNaP = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_NaP = dmNaP(dyna.m_NaP, p.dt*1e3) + dyna.m_NaP
        dyna.h_NaP = dhNaP(dyna.h_NaP, p.dt*1e3) + dyna.h_NaP

        # calculate the open-probability of the channel:
        P = (dyna.m_NaP ** self._mpower) * (dyna.h_NaP ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg2[sim.iNa][dyna.targets_vgNaP] = dyna.maxDmNaP * P


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

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the persistent vgNa channel: Nav1p6')

        self._mpower = 1.0
        self._hpower = 0.0

        dyna.m_NaP = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        dyna.h_NaP = 1.0

    def _calculate_state(self, V, dyna, sim, p):

        self._mInf = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self._mTau = 1

        self._hInf = 1.0
        self._hTau = 1.0