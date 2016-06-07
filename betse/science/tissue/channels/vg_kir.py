#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated inward rectifying potassium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np

from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb

# FIXME update with new formalism dealing with charge transfer


# .................... BASE                               ....................
class VgKABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated inward rectifying potassium channel classes.

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
        Initialize targeted voltage-gated inward rectifying channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        V = sim.vm[dyna.targets_vgKir] * 1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        self._calculate_state(
            V=sim.vm[dyna.targets_vgKir] * 1000,
            dyna=dyna, sim=sim, p=p)

        self._implement_state(dyna, sim, cells, p)

    def _implement_state(self, dyna, sim, cells, p):
        # calculate m and h channel states using RK4:
        dmKir = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhKir = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Kir = dmKir(dyna.m_Kir, p.dt * 1e3) + dyna.m_Kir
        dyna.h_Kir = dhKir(dyna.h_Kir, p.dt * 1e3) + dyna.h_Kir

        # calculate the open-probability of the channel:
        P = (dyna.m_Kir ** self._mpower) * (dyna.h_Kir ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg2[sim.iK][dyna.targets_vgKir] = dyna.maxDmKir * P


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
class Kir2p1(VgKABC):  # FIXME do we have to prefix these by vgK -- how about just Kir2p1?
    '''
    Kir 2.1 model from Makary et al.

    This channel has a greater tendency to allow potassium to flow into a cell rather than out of a cell,
    probably participates in establishing action potential waveform and excitability of neuronal and muscle tissues.
    Mutations in this gene have been associated with Andersen syndrome, which is characterized by periodic paralysis,
    cardiac arrhythmias, and dysmorphic features

    Reference: Makary SM. et al. A difference in inward rectification and polyamine block and permeation between
     the Kir2.1 and Kir3.1/Kir3.4 K+ channels. J. Physiol. (Lond.), 2005 Nov 1 , 568 (749-66).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the inward rectifying K+ channel: Kir2p1')

        self.vrev = -70.6     # reversal voltage used in model [mV]
        Texpt = 28    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0   # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Kir = 1 / (1 + np.exp((V - (-96.48)) / 23.26))
        dyna.h_Kir = 1 / (1 + np.exp((V - (-168.28)) / -44.13))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 2


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - (-96.48)) / 23.26))
        self._mTau = 3.7 + (-3.37 / (1 + np.exp((V - -32.9) / 27.93)))
        self._hInf = 1 / (1 + np.exp((V - (-168.28)) / -44.13))
        self._hTau = 0.85 + (306.3 / (1 + np.exp((V - -118.29) / -27.23)))
