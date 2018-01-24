#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated calcium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channelo.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

# .................... BASE                               ....................
class VgCaABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated calcium channel classes.

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
        Initialize targeted voltage-gated calcium channels at the initial
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

        self.ions = ['Ca']
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
        Do something.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V):
        '''
        Do something.
        '''
        pass

# ....................{ SUBCLASS                           }....................

class Cav3p3(VgCaABC):
    '''

    T-type calcium channel model from Traboulsie et al. Based on a Cav3.3 channel.


    Reference: Traboulsie A. et al. Subunit-specific modulation of T-type calcium
    channels by zinc. J. Physiol. (Lond.), 2007 Jan 1 , 578 (159-71).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Cav3p3')

        self.time_unit = 1.0e3

        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp((V - -45.454426) / -5.073015))
        self.h = 1 / (1 + np.exp((V - (-74.031965)) / 8.416382))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - -45.454426) / -5.073015))
        self._mTau = 3.394938 + (54.187616 / (1 + np.exp((V - -40.040397) / 4.110392)))
        self._hInf = 1 / (1 + np.exp((V - (-74.031965)) / 8.416382))
        self._hTau = 109.701136 + (0.003816 * np.exp(-V / 4.781719))

class Cav2p1(VgCaABC):
    '''

    P/Q-type calcium channel (Cav2.1) model from Miyasho et al.

    Found in the nervous system, especially Purkunjie neurons,  play an important role in mediating neurotransmitter
    release in the nervous system, postsynaptic integration, neuroplasticity, neural excitability,
    and gene transcription.


    Reference: Miyasho T. et al. Low-threshold potassium channels and a low-threshold calcium channel
    regulate Ca2+ spike firing in the dendrites of cerebellar Purkinje neurons: a modeling study.
    Brain Res., 2001 Feb 9 , 891 (106-15).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Cav2p1')

        self.time_unit = 1.0e3

        self.vrev = 135.0     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        mAlpha = 8.5 / (1 + np.exp((V - 8) / (-12.5)))
        mBeta = 35 / (1 + np.exp((V + 74) / (14.5)))

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = mAlpha / (mAlpha + mBeta)
        self.h = 1.0

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        mAlpha = 8.5 / (1 + np.exp((V - 8) / (-12.5)))
        mBeta = 35 / (1 + np.exp((V + 74) / (14.5)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        self._hInf = 1.0
        self._hTau = 1.0

class Cav1p3(VgCaABC):
    '''

    Low voltage activating L-type calcium channel model from Avery et al.


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Cav1.3')

        self.time_unit = 1.0e3

        self.vrev = 113.0     # reversal voltage used in model [mV]
        Texpt = 21.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)

        # V = V - 10

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self.h = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = (5.0000 + 20.0000 / (1 + np.exp((V - -25.000) / 5)))
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = (20.0000 + 50.0000 / (1 + np.exp((V - -40.000) / 7)))

class Cav1p2(VgCaABC):
    '''

    L-type calcium channel model of Carlin et. al.

    Reference: Carlin KP. et al. Characterization of calcium
    currents in functionally mature mouse spinal motoneurons. Eur. J. Neurosci., 2000 May , 12 (1624-34).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Cav1.2')

        self.time_unit = 1.0e3

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 21.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)

        V = V - 10

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self.h = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        V = V - 10
        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = (5.0000 + 20.0000 / (1 + np.exp((V - -25.000) / 5)))
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = (20.0000 + 50.0000 / (1 + np.exp((V - -40.000) / 7)))

class Cav2p3(VgCaABC):
    '''

    R-type calcium channel model from Miyasho et al. based on Cav2p3.


    Reference: Miyasho T. et al. Low-threshold potassium channels and a low-threshold calcium channel
    regulate Ca2+ spike firing in the dendrites of cerebellar Purkinje neurons:
    a modeling study. Brain Res., 2001 Feb 9 , 891 (106-15).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the R-type vgCa channel type: Cav2p3')

        self.time_unit = 1.0e3

        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        mAlpha = 2.6 / (1 + np.exp((V + 7) / -8))
        mBeta = 0.18 / (1 + np.exp((V + 26) / 4))
        hAlpha = 0.0025 / (1 + np.exp((V + 32) / 8))
        hBeta = 0.19 / (1 + np.exp((V + 42) / -10))

        self.m = mAlpha / (mAlpha + mBeta)
        self.h = hAlpha / (hAlpha + hBeta)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        mAlpha = 2.6 / (1 + np.exp((V + 7) / -8))
        mBeta = 0.18 / (1 + np.exp((V + 26) / 4))
        hAlpha = 0.0025 / (1 + np.exp((V + 32) / 8))
        hBeta = 0.19 / (1 + np.exp((V + 42) / -10))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf =hAlpha / (hAlpha + hBeta)
        self._hTau = 1 / (hAlpha + hBeta)

class Cav2p2(VgCaABC):
    '''

    N-type calcium channel model from Huang SJ et al.


    Reference: Huang SJ. et al. Activation and inactivation properties of voltage-gated calcium currents
    in developing cat retinal ganglion cells. Neuroscience, 1998 Jul , 85 (239-47).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the N-type vgCa channel type: Cav2p2')

        self.time_unit = 1.0e3

        self.vrev = 135     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0

        # indsV = (V == 20.0).nonzero()
        # V[indsV] += 1.0e-6

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        mAlpha = (0.1 * (V - 20) / (1 - np.exp(-(V - 20) / 10)))
        mBeta = 0.4 * np.exp(-(V + 25) / 18)
        hAlpha = 0.01 * np.exp(-(V + 50) / 10)
        hBeta =  0.1 / (1 + np.exp(-(V + 17) / 17))


        self.m = mAlpha / (mAlpha + mBeta)
        self.h = hAlpha / (hAlpha + hBeta)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # indsV = (V == 20.0).nonzero()
        # V[indsV] += 1.0e-6

        mAlpha = (0.1 * (V - 20) / (1 - np.exp(-(V - 20) / 10)))
        mBeta = 0.4 * np.exp(-(V + 25) / 18)
        hAlpha = 0.01 * np.exp(-(V + 50) / 10)
        hBeta =  0.1 / (1 + np.exp(-(V + 17) / 17))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf =hAlpha / (hAlpha + hBeta)
        self._hTau = 1 / (hAlpha + hBeta)

class Cav3p1(VgCaABC):
    '''

    T-type calcium channel model from Traboulsie et al, based on Cav3p1.


    Reference: Traboulsie A. et al. Subunit-specific modulation of T-type calcium channels by zinc.
    J. Physiol. (Lond.), 2007 Jan 1 , 578 (159-71).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the T-type vgCa channel type: Cav3p1')

        self.time_unit = 1.0e3

        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp((V - (-42.921064)) / -5.163208))
        self.h = 1 / (1 + np.exp((V - (-72.907420)) / 4.575763))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - (-42.921064)) / -5.163208))

        inds_lt10 = (V < -10e-3).nonzero()
        inds_gt10 = (V >= -10e-3).nonzero()

        self._mTau = -0.855809 + (1.493527 * np.exp(-V / 27.414182))
        self._mTau[inds_gt10] = 1.0
        self._hInf = 1 / (1 + np.exp((V - (-72.907420)) / 4.575763))
        self._hTau = 9.987873 + (0.002883 * np.exp(-V / 5.598574))

class Ca_L2(VgCaABC):
    '''

    L-type calcium channel model Avery et al. L-type channels are higher-voltage activating and very
    persistent. They are commonly found in muscle or glands, where they induce activities such as hormone
    release or muscle contraction in response to neural stimulation.

    This channel has been modified to activate at more depolarized Vmem.


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L2')

        self.time_unit = 1.0e3

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self.h = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = 10.0
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = 59.0

class Ca_L3(VgCaABC):
    '''

    L-type calcium channel model Avery et al. L-type channels are higher-voltage activating and very
    persistent. They are commonly found in muscle or glands, where they induce activities such as hormone
    release or muscle contraction in response to neural stimulation.

    This channel has been modified to activate at more depolarized Vmem and to have a broader activation peak.


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L3')

        self.time_unit = 1.0e3

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0

        V = V - 15

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self.h = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        # self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        # self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        # self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))

        V = V - 15

        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = 10.0
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = 59.0

class Cav_G(VgCaABC):
    '''

    General calcium channel model from Reuveni et al.

    Reference: Reuveni I. et al. Stepwise repolarization from Ca2+ plateaus in neocortical pyramidal cells: evidence
    for nonhomogeneous distribution of HVA Ca2+ channels in dendrites. J. Neurosci., 1993 Nov , 13 (4609-21).




    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Cav_G')

        self.time_unit = 1.0e3

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0


        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        mAlpha = (0.055 * (-27 - V)) / (np.exp((-27 - V) / 3.8) - 1)
        mBeta = (0.94 * np.exp((-75 - V) / 17))
        self.m = mAlpha / (mAlpha + mBeta)

        hAlpha = (0.000457 * np.exp((-13 - V) / 50))
        hBeta = (0.0065 / (np.exp((-V - 15) / 28) + 1))
        self.h = hAlpha / (hAlpha + hBeta)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """


        mAlpha = (0.055 * (-27 - V)) / (np.exp((-27 - V) / 3.8) - 1)
        mBeta = (0.94 * np.exp((-75 - V) / 17))
        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        hAlpha = (0.000457 * np.exp((-13 - V) / 50))
        hBeta = (0.0065 / (np.exp((-V - 15) / 28) + 1))
        self._hInf = hAlpha / (hAlpha + hBeta)
        self._hTau = 1 / (hAlpha + hBeta)

class CaLeak(VgCaABC):

    '''
    Simple calcium leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a substance-modulated Ca++ channel')

        self.time_unit = 1.0

        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        self.qt = 1.0

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




