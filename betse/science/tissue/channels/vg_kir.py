#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
from betse.science import sim_toolbox as stb

# FIXME update with new formalism dealing with charge transfer


# .................... BASE                               ....................
class VgKirABC(ChannelsABC, metaclass=ABCMeta):
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

        self.modulator = 1.0

        self.v_corr = 0.0

        V = sim.vm*1000 + self.v_corr

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        V = sim.vm*1000 + self.v_corr

        self._calculate_state(
            V=V,
            dyna=dyna, sim=sim, p=p)

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):

        # calculate m and h channel states using RK4:
        dmKir = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhKir = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Kir = dmKir(dyna.m_Kir, p.dt * 1e3) + dyna.m_Kir
        dyna.h_Kir = dhKir(dyna.h_Kir, p.dt * 1e3) + dyna.h_Kir

        # calculate the open-probability of the channel:
        P = (dyna.m_Kir ** self._mpower) * (dyna.h_Kir ** self._hpower)

        # print(P.min(), P.max(), P.mean())

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem = sim.cc_cells[sim.iK][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env = sim.cc_env[sim.iK][cells.map_mem2ecm]

        else:
            c_env = sim.cc_env[sim.iK]

        IdM = np.ones(sim.mdl)

        z_ion = sim.zs[sim.iK] * IdM

        # make sure probablility is binned between zero and 1:
        # inds_PL = (P < 0.0).nonzero()
        # inds_PH = (P > 1.0).nonzero()
        #
        # P[inds_PL] = 0.0
        # P[inds_PH] = 1.0

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmKir*P*1.0e-9

        # calculate specific ion flux contribution for this channel:
        delta_Q = stb.electroflux(c_env, c_mem, Dchan, p.tm * IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iK, delta_Q, dyna.targets_vgKir, sim, cells, p)


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
class Kir2p1(VgKirABC):
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

        # self.rectification = -1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - (-96.48)) / 23.26))
        self._mTau = 3.7 + (-3.37 / (1 + np.exp((V - -32.9) / 27.93)))
        self._hInf = 1 / (1 + np.exp((V - (-168.28)) / -44.13))
        self._hTau = 0.85 + (306.3 / (1 + np.exp((V - -118.29) / -27.23)))

class Kv3p4(VgKirABC):

    """

     Kv3p4 is a rapidly inactivating, A-type inward/outward potassium current.

     Reference: Vega-Saenz de Miera E. et al. Cloning of ShIII (Shaw-like) cDNAs encoding a novel
     high-voltage-activating, TEA-sensitive, type-A K+ channel. Proc. Biol. Sci., 1992 Apr 22 , 248 (9-18).

    """

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the A-current K+ channel: Kv3p4')

        self.vrev = -65.0  # reversal voltage used in model [mV]
        Texpt = 28  # temperature of the model in degrees C
        simT = sim.T - 273  # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Kir = 1 / (1 + np.exp(((V - (-3.400)) / (-8.400))))
        dyna.h_Kir =  1 / (1 + np.exp(((V - (-53.320)) / (7.400))))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp(((V - (-3.400)) / (-8.400))))
        self._mTau = 10.000 / (1 + np.exp(((V - (4.440)) / (38.140))))
        self._hInf = 1 / (1 + np.exp(((V - (-53.320)) / (7.400))))
        self._hTau = 20000.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))

class Kv3p3(VgKirABC):


    """

     Kv3p3 is a rapidly inactivating, A-type inward/outward potassium current, implicated in rapid firing of
     neurons and is primarily found in Purkunjie neurons.

     Reference: Desai R. et al. Protein kinase C modulates inactivation of Kv3.3 channels.
     J. Biol. Chem., 2008 Aug 8 , 283 (22283-94).

    """

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the A-current K+ channel: Kv3p3')

        self.vrev = 82.0  # reversal voltage used in model [mV]
        Texpt = 28  # temperature of the model in degrees C
        simT = sim.T - 273  # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Kir = 1 / (1 + np.exp((V - 35) / -7.3))
        dyna.h_Kir =  0.25 + (0.75 / (1 + np.exp((V - (-28.293856)) / 29.385636)))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):

        self._mInf = 1 / (1 + np.exp((V - 35) / -7.3))
        self._mTau = 0.676808 + (27.913114 / (1 + np.exp((V - 22.414149) / 9.704638)))
        self._hInf = 0.25 + (0.75 / (1 + np.exp((V - (-28.293856)) / 29.385636)))
        self._hTau = 199.786728 + (2776.119438 * np.exp(-V / 7.309565))


