#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated calcium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np
from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb


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

    def init(self, dyna, sim, cells, p):
        '''
        Initialize targeted voltage-gated calcium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        self.modulator = 1.0 # init the modulator to 1.0, which may be updated programmatically later

        self.v_corr = 0.0  # offset of voltages in the model -- experimental junction voltage [mV]

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

        self._calculate_state(V, dyna=dyna, sim=sim, p=p)

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):

        # calculate m and h channel states using RK4:
        dmCa = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhCa = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Ca = dmCa(dyna.m_Ca, p.dt*1e3) + dyna.m_Ca
        dyna.h_Ca = dhCa(dyna.h_Ca, p.dt*1e3) + dyna.h_Ca

        # calculate the open-probability of the channel:
        P = (dyna.m_Ca ** self._mpower) * (dyna.h_Ca ** self._hpower)



        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the BETSE electrodiffusion equation (GHK flux equation).

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        # delta_Q = -(dyna.maxDmCa*P*(V - self.vrev))


        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem = sim.cc_cells[sim.iCa][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env = sim.cc_env[sim.iCa][cells.map_mem2ecm]

        else:
            c_env = sim.cc_env[sim.iCa]

        IdM = np.ones(sim.mdl)

        z_ion = sim.zs[sim.iCa] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmCa*P*1.0e-9

        # calculate specific ion flux contribution for this channel:
        delta_Q = stb.electroflux(c_env, c_mem, Dchan, p.tm * IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

        # the cube power in the vgNa expression is rather difficult mathematically, but necessary
        # clip the unreasonably high portions of the Na+ flux, so as not to overload the system:
        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iCa, delta_Q, dyna.targets_vgCa, sim, cells, p)



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
class Ca_L(VgCaABC):
    '''

    L-type calcium channel model Avery et al. L-type channels are higher-voltage activating and very
    persistent. They are commonly found in muscle or glands, where they induce activities such as hormone
    release or muscle contraction in response to neural stimulation.


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L')

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        dyna.h_Ca = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        # self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        # self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        # self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))

        # V = V + 10
        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = (5.0000 + 20.0000 / (1 + np.exp((V - -25.000) / 5)))
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = (20.0000 + 50.0000 / (1 + np.exp((V - -40.000) / 7)))
        # V = V - 10

class Ca_T(VgCaABC):
    '''

    T-type calcium channel model from Traboulsie et al. Based on a Cav3.3 channel.


    Reference: Traboulsie A. et al. Subunit-specific modulation of T-type calcium
    channels by zinc. J. Physiol. (Lond.), 2007 Jan 1 , 578 (159-71).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_T')

        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1 / (1 + np.exp((V - -45.454426) / -5.073015))
        dyna.h_Ca = 1 / (1 + np.exp((V - (-74.031965)) / 8.416382))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - -45.454426) / -5.073015))
        self._mTau = 3.394938 + (54.187616 / (1 + np.exp((V - -40.040397) / 4.110392)))
        self._hInf = 1 / (1 + np.exp((V - (-74.031965)) / 8.416382))
        self._hTau = 109.701136 + (0.003816 * np.exp(-V / 4.781719))

class Ca_PQ(VgCaABC):
    '''

    P/Q-type calcium channel (Cav2.1) model from Miyasho et al.

    Found in the nervous system, especially Purkunjie neurons,  play an important role in mediating neurotransmitter
    release in the nervous system, postsynaptic integration, neuroplasticity, neural excitability,
    and gene transcription.


    Reference: Miyasho T. et al. Low-threshold potassium channels and a low-threshold calcium channel
    regulate Ca2+ spike firing in the dendrites of cerebellar Purkinje neurons: a modeling study.
    Brain Res., 2001 Feb 9 , 891 (106-15).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_PQ')

        self.vrev = 135.0     # reversal voltage used in model [mV]
        Texpt = 0.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        mAlpha = 8.5 / (1 + np.exp((V - 8) / (-12.5)))
        mBeta = 35 / (1 + np.exp((V + 74) / (14.5)))

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = mAlpha / (mAlpha + mBeta)
        dyna.h_Ca = 1.0

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
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

class Ca_L2(VgCaABC):
    '''

    L-type calcium channel model Carlin et al. L-type channels are higher-voltage activating and very
    persistent. They are commonly found in muscle or glands, where they induce activities such as hormone
    release or muscle contraction in response to neural stimulation.


    Reference: Carlin KP. et al. Characterization of calcium currents in functionally mature mouse spinal
    motoneurons. Eur. J. Neurosci., 2000 May , 12 (1624-34).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L2')

        self.vrev = 80.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        dyna.h_Ca = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        # self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        # self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        # self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))
        self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))


class Ca_L_lv(VgCaABC):
    '''

    Low voltage activating L-type calcium channel model altered from Avery et al.


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L_lv')

        self.vrev = 113.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # V = V + 50

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        dyna.h_Ca = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        # self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        # self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        # self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))

        V = V + 10
        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = (5.0000 + 20.0000 / (1 + np.exp((V - -25.000) / 5)))
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = (20.0000 + 50.0000 / (1 + np.exp((V - -40.000) / 7)))

class Ca_L3(VgCaABC):
    '''

    L-type calcium channel model Avery et al. L-type channels are higher-voltage activating and very
    persistent. They are commonly found in muscle or glands, where they induce activities such as hormone
    release or muscle contraction in response to neural stimulation.

    This channel has been modified to activate at higher Vmem (peak at about -10 mV)


    Reference: Avery RB. et al. Multiple channel types contribute to the low-voltage-activated calcium current in
    hippocampal CA3 pyramidal neurons. J. Neurosci., 1996 Sep 15 , 16 (5567-82).

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgCa channel type: Ca_L3')

        self.vrev = 131.0     # reversal voltage used in model [mV]
        Texpt = 36.0    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        V = V - 20

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        dyna.h_Ca = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        # self._mInf = 1.0000 / (1 + np.exp(((V - 10) + 30.000) / -6))
        # self._mTau = 5.0000 + 20.0000 / (1 + np.exp(((V - 10) + 25.000) / 5))
        # self._hInf = 1.0000 / (1 + np.exp(((V - 10) + 80.000) / 6.4))
        # self._hTau = 20.0000 + 50.0000 / (1 + np.exp(((V - 10) + 40.000) / 7))

        V = V - 20
        self._mInf = 1.0000 / (1 + np.exp((V - -30.000) / -6))
        self._mTau = (5.0000 + 20.0000 / (1 + np.exp((V - -25.000) / 5)))
        self._hInf = 1.0000 / (1 + np.exp((V - -80.000) / 6.4))
        self._hTau = (20.0000 + 50.0000 / (1 + np.exp((V - -40.000) / 7)))
        # V = V - 10


class CaLeak(VgCaABC):

    '''
    Simple calcium leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a substance-modulated Ca++ channel')


        self.vrev = 30     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        self.qt = 1.0

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Ca = 1
        dyna.h_Ca = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1
        self._mTau = 1
        self._hInf = 1
        self._hTau = 1




