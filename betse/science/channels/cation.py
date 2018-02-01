#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Cation leak channel classes.
'''

# .................... IMPORTS                            ....................
import numpy as np
from abc import ABCMeta, abstractmethod
from betse.science import sim_toolbox as stb
from betse.science.channels.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs


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

    def init(self, vm, cells, p, targets = None):
        '''
        Initialize targeted voltage-gated funny current channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
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

        self.ions = ['Na', 'K', 'Ca']
        self.rel_perm = [self._PmNa, 1.0, self._PmCa]


    def run(self, vm, p):
        '''
        Handle all targeted voltage-gated funny current channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        V = vm[self.targets]*1000

        self._calculate_state(V)

        self._implement_state(V, p)

    def _implement_state(self, V, p):

        # Update the 'm' and 'h' gating functions:
        self.update_mh(p, time_unit = self.time_unit)

        # calculate the open-probability of the channel:
        P = (self.m ** self._mpower) * (self.h ** self._hpower)

        self.P = np.zeros(self.mdl)
        self.P[self.targets] = P


    @abstractmethod
    def _init_state(self, V):
        '''
        Initializes values of the m and h gates of the channel.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V):
        '''
        Calculates time-dependent values of the m and h gates of the channel.
        '''
        pass

# ....................{ SUBCLASS                           }....................
class CatLeak(CationABC):
    '''

    Membrane leak channel letting Na+, K+ and Ca2+ into the cell

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the cation leak channel')

        self.time_unit = 1.0

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = np.ones(self.data_length)
        self.h = np.ones(self.data_length)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0

        # ratio of various ion membrane permeability of channel (relative to K+)
        self._pmNa = 1.0
        self._pmCa = 0.0


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0


class CatLeak2(CationABC):
    '''

    Membrane leak channel letting Na+, K+ and Ca2+ into the cell

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the cation leak channel Na 2:K 1')

        self.time_unit = 1.0

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = np.ones(self.data_length)
        self.h = np.ones(self.data_length)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0

        # ratio of Na to K membrane permeability of channel
        self._pmNa = 2.0
        self._pmCa = 0.0


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0