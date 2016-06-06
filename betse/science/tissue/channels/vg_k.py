#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated potassium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np

from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb


# .................... BASE                               ....................
class VgKABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated potassium channel classes.

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

        V = sim.vm[dyna.targets_vgK] * 1000

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
            V=sim.vm[dyna.targets_vgK] * 1000,
            dyna=dyna, sim=sim, p=p)

        self._implement_state(dyna, sim, p)

    def _implement_state(self, dyna, sim, p):
        # calculate m and h channel states using RK4:
        dmK = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhK = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_K = dmK(dyna.m_K, p.dt * 1e3) + dyna.m_K
        dyna.h_K = dhK(dyna.h_K, p.dt * 1e3) + dyna.h_K

        # calculate the open-probability of the channel:
        P = (dyna.m_K ** self._mpower) * (dyna.h_K ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iK][dyna.targets_vgK] = dyna.maxDmK * P


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
class vgK_Default(VgKABC):
    '''
    Kv1.2 model from  Sprunger et al.

    This channel produces well behaved action-potentials with a variety of vgNa channels. Good general-purpose
    vgK channel.

    Potassium voltage-gated channel Kv1.2 is a member of the shaker-related subfamily and belongs to the
    delayed rectifier class of channels, which allow nerve cells to efficiently re-polarize following an
    action potential. Kv1.2 is primarily found in the central nervous system.

    Reference: Sprunger LK. et al. Effects of charybdotoxin on K+ channel (KV1.2)
    deactivation and inactivation kinetics. Eur. J. Pharmacol., 1996 Oct 31 , 314 (357-64).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel type: vgK_Default (Kv1.2)')

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_K = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
        dyna.h_K = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
        self._mTau = 150.0000/(1+ np.exp((V + 67.5600)/34.1479))
        self._hInf = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))
        self._hTau = 15000.0000/(1+ np.exp(-(V + 46.5600)/44.1479))

class vgK_Kv1p1(VgKABC):
    '''
    Kv1.1 model from Christie et al, cloned from rat brain, studied in xenopus.

    Kv1.1 are low-voltage activated (LVA) channels, expressed primarily in the central nervous system, and
    which open with small depolarizations at or below resting potential

    Reference: Christie MJ. et al. Expression of a cloned rat brain potassium channel in Xenopus oocytes.
    Science, 1989 Apr 14 , 244 (221-4).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel type: vgK_Kv1p1 ')

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_K = 1.0000 / (1 + np.exp((V - -30.5000) / -11.3943))
        dyna.h_K = 1.0000 / (1 + np.exp((V - -30.0000) / 27.3943))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 2


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -30.5000) / -11.3943))
        self._mTau = 30.0000 / (1 + np.exp((V - -76.5600) / 26.1479))
        self._hInf = 1.0000 / (1 + np.exp((V - -30.0000) / 27.3943))
        self._hTau = 15000.0000 / (1 + np.exp((V - -160.5600) / -100.0000))

class vgK_Kv1p5(VgKABC):
    '''
    Kv1.5 model from

    This channel produces well behaved action-potentials with a variety of vgNa channels. Good general-purpose
    vgK channel.

    Kv1.5 underlies the cardiac ultra-rapid delayed rectifier potassium current (IKur) in humans.
    In addition, Kv1.5 channels are also expressed in many other organs, such as pulmonary arteries,
    brain, skeletal muscle, and have crucial function in cell cycle regulation.

    Reference: Philipson LH. et al. Sequence and functional expression in Xenopus oocytes of a
    human insulinoma and islet potassium channel. Proc. Natl. Acad. Sci. U.S.A., 1991 Jan 1 , 88 (53-7).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel type: vgK_Kv1p5 ')

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        dyna.m_K = 1.0000 / (1 + np.exp((V - -6.0000) / -6.4000))
        dyna.h_K = 1.0000 / (1 + np.exp((V - -25.3000) / 3.5000))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -6.0000) / -6.4000))
        self._mTau = (-0.1163 * V) + 8.3300
        self._hInf = 1.0000 / (1 + np.exp((V - -25.3000) / 3.5000))
        self._hTau = (-15.5000 * V) + 1620.0000


