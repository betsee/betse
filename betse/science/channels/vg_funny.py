#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated funny current channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channels.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

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

    def init(self, vm, cells, p, targets = None):
        '''
        Initialize targeted voltage-gated funny current channels at the initial
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

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN2')

        self.time_unit = 1.0e3

        self.v_corr = 0

        V = V - 10

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -99) / 6.2))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability to Ca2+
        self._PmNa = 0.2  # channel permeability ratio to Na+


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        V = V - 10

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

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN4')

        self.time_unit = 1.0e3

        self.v_corr = 0

        V = V - 10

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -100) / 9.6))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability to Ca2+
        self._PmNa = 0.2  # channel permeability ratio to Na+


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        V = V - 10

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

    def _init_state(self, V):

        logs.log_info('You are using the funny current channel: HCN1')

        self.time_unit = 1.0e3

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -94) / 8.1))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability ratio to Ca2+
        self._PmNa = 0.2  # channel permeability ratio to Na+

    def _calculate_state(self, V):


        self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0000 / (1 + np.exp((V - -94) / 8.1))
        self._mTau = 30.0000
        self._hInf = 1
        self._hTau = 1

class HCNLeak(VgFunABC):
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

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN Leak')

        self.time_unit = 1.0

        self.v_corr = 0

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = np.ones(self.data_length)
        self.h = np.ones(self.data_length)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability to Ca2+
        self._PmNa = 0.33  # channel permeability ratio to Na+


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0

class HCN2_cAMP(VgFunABC):
    '''
    HCN2 model from 21 day old dorsal root ganglion of mouse. HCN channels are voltage-gated ionic channels,
    regulated by cyclic nucleotides, such as cyclic adenosine-mono-phosphate (cAMP), modelled here as a simple
    shift in the channel voltage gating threshold.
    In contrast to most Na+ and K+ ionic channels, which open when membrane potential is depolarized,
    they are opened when the membrane potential hyperpolarizes below -50 mV.

    This channel is expressed in brain and heart tissue, though the specific function is unknown.

    Reference:  Moosmang S. et al. Cellular expression and functional characterization of four
    hyperpolarization-activated pacemaker channels in cardiac and neuronal tissues.
    Eur. J. Biochem., 2001 Mar , 268 (1646-52).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN2 with cAMP activation.')

        self.time_unit = 1.0e3

        self.v_corr = 0

        V = V - 20

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V + 99) / 6.2))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability to Ca2+
        self._PmNa = 0.2  # channel permeability ratio to Na+


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        V = V - 20

        self._mInf = 1.0000 / (1 + np.exp((V + 99) / 6.2))
        self._mTau = 184.0000
        self._hInf = 1
        self._hTau = 1

class HCN4_cAMP(VgFunABC):
    '''
    HCN4 model from 21 day old dorsal root ganglion of mouse. HCN channels are voltage-gated ionic channels,
    regulated by cyclic nucleotides, such as cyclic adenosine-mono-phosphate (cAMP), modelled here as a simple
    shift in the channel voltage gating threshold.
    In contrast to most Na+ and K+ ionic channels, which open when membrane potential is depolarized,
    they are opened when the membrane potential hyperpolarizes below -50 mV.

    This channel is expressed in brain and heart tissue, and plays a crucial role in pacemaker activity of the
    heart, and possibly the repetitive firing of Purkunjie neurons.

    Reference:  Moosmang S. et al. Cellular expression and functional characterization of four
    hyperpolarization-activated pacemaker channels in cardiac and neuronal tissues.
    Eur. J. Biochem., 2001 Mar , 268 (1646-52).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the funny current channel: HCN4 with cAMP activation.')

        self.time_unit = 1.0e3

        self.v_corr = 0

        V = V - 20

        # initialize values of the m and h gates of the HCN2 based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V  + 100) / 9.6))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

        self._PmCa = 0.05  # channel permeability to Ca2+
        self._PmNa = 0.20  # channel permeability ratio to Na+


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = -45  # reversal voltage used in model [mV]

        V = V - 24

        self._mInf = 1.0000 / (1 + np.exp((V + 100) / 9.6))
        self._mTau = 461.0000
        self._hInf = 1
        self._hTau = 1

