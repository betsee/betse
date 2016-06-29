#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated funny current channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np

from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb


# .................... BASE                               ....................
class VgFunABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated funny current channel classes.

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

        V = sim.vm[dyna.targets_vgFun] * 1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated funny current channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        V = sim.vm[dyna.targets_vgFun] * 1000

        self._calculate_state(V, dyna, sim, p)

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):
        # calculate m and h channel states using RK4:
        dmFun = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhFun = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Fun = dmFun(dyna.m_Fun, p.dt * 1e3) + dyna.m_Fun
        dyna.h_Fun = dhFun(dyna.h_Fun, p.dt * 1e3) + dyna.h_Fun

        # calculate the open-probability of the channel:
        P = (dyna.m_Fun ** self._mpower) * (dyna.h_Fun ** self._hpower)

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        delta_Q = - (dyna.maxDmFun*P*(V - self.vrev))

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q, dyna.targets_vgFun, sim, cells, p)
        self.update_charge(sim.iK, delta_Q, dyna.targets_vgFun, sim, cells, p)


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
class HCN2(VgFunABC):
    '''
    HCN2 model from 21 day old dorsal root ganglion of mouse. HCN channels are voltage-gated ionic channels,
    regulated by cyclic nucleotides, such as cyclic adenosine-mono-phosphate (cAMP) (not modelled here).
    In contrast to most Na+ and K+ ionic channels, which open when membrane potential is depolarized,
    they are opened when the membrane potential hyperpolarizes below -50 mV.

    This channel is expressed in brain and heart tissue, though the specific function is unknown.

    Reference:  Moosmang S. et al. Cellular expression and functional characterization of four
    hyperpolarization-activated pacemaker channels in cardiac and neuronal tissues.
    Eur. J. Biochem., 2001 Mar , 268 (1646-52).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN2')

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        dyna.m_Fun = 1.0000 / (1 + np.exp((V - -99) / 6.2))
        dyna.h_Fun = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0000 / (1 + np.exp((V - -99) / 6.2))
        self._mTau = 184.0000
        self._hInf = 1
        self._hTau = 1

class HCN4(VgFunABC):
    '''
    HCN4 model from 21 day old dorsal root ganglion of mouse. HCN channels are voltage-gated ionic channels,
    regulated by cyclic nucleotides, such as cyclic adenosine-mono-phosphate (cAMP) (not modelled here).
    In contrast to most Na+ and K+ ionic channels, which open when membrane potential is depolarized,
    they are opened when the membrane potential hyperpolarizes below -50 mV.

    This channel is expressed in brain and heart tissue, and plays a crucial role in pacemaker activity of the
    heart, and possibly the repetitive firing of Purkunjie neurons.

    Reference:  Moosmang S. et al. Cellular expression and functional characterization of four
    hyperpolarization-activated pacemaker channels in cardiac and neuronal tissues.
    Eur. J. Biochem., 2001 Mar , 268 (1646-52).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN4')

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        dyna.m_Fun = 1.0000 / (1 + np.exp((V - -100) / 9.6))
        dyna.h_Fun = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0000 / (1 + np.exp((V - -100) / 9.6))
        self._mTau = 461.0000
        self._hInf = 1
        self._hTau = 1

class HCN1(VgFunABC):

    """
    Fastest activation of the HCN channels, the HCN1 channel is also implicated in sinus rhythm control
    of the heart, brain activity, and in the retina.

    Reference: 3] 	Moosmang S. et al. Cellular expression and functional characterization of four
    hyperpolarization-activated pacemaker channels in cardiac and neuronal tissues.
    Eur. J. Biochem., 2001 Mar , 268 (1646-52).
    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the funny current channel: HCN1')

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        dyna.m_Fun = 1.0000 / (1 + np.exp((V - -94) / 8.1))
        dyna.h_Fun = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

    def _calculate_state(self, V, dyna, sim, p):


        self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0000 / (1 + np.exp((V - -94) / 8.1))
        self._mTau = 30.0000
        self._hInf = 1
        self._hTau = 1