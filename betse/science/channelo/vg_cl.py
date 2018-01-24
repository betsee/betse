#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated chloride channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channelo.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

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

    def init(self, vm, cells, p, targets = None):
        '''
        Initialize targeted voltage-gated potassium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        if targets is None:

            self.targets = None

            if type(vm) == float or type(vm) == np.float64:
                self.data_length = 1

            else:
                self.data_length = len(vm)


        else:
            self.targets = targets
            self.data_length = len(self.targets)
            self.mdl = len(cells.mem_i)

        self.modulator = 1.0

        self.v_corr = 0.0  # in experiments, the measurement junction voltage is about 10 mV

        if self.targets is None:

            V = vm * 1000 + self.v_corr

        else:
            V = vm[self.targets] * 1000 + self.v_corr


        self._init_state(V)

        self.ions = ['Cl']
        self.rel_perm = [1.0]


    def run(self, vm, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        if self.targets is None:

            V = vm*1000 + self.v_corr

        else:

            V = vm[self.targets] * 1000 + self.v_corr

        self._calculate_state(V)

        self._implement_state(V, p)

    def _implement_state(self, V, p):

        # Update the 'm' and 'h' gating functions:
        self.update_mh(p, time_unit = self.time_unit)

        # calculate the open-probability of the channel:
        P = (self.m ** self._mpower) * (self.h ** self._hpower)

        if self.targets is None:

            self.P = P

        else:
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


class ClLeak(VgClABC):

    '''
    Simple chloride leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a Cl- leak channel')

        self.time_unit = 1.0


        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = np.ones(self.data_length)
        self.h = np.ones(self.data_length)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1
        self._mTau = 1
        self._hInf = 1
        self._hTau = 1

