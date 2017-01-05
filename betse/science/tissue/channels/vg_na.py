#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
from betse.science import sim_toolbox as stb


# .................... BASE                               ....................
class VgNaABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated sodium channel classes.

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
        Initialize targeted voltage-gated sodium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        self.modulator = 1.0

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
        dmNa = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhNa = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Na = dmNa(dyna.m_Na, p.dt*1e3) + dyna.m_Na
        dyna.h_Na = dhNa(dyna.h_Na, p.dt*1e3) + dyna.h_Na

        # calculate the open-probability of the channel:
        P = (dyna.m_Na ** self._mpower) * (dyna.h_Na ** self._hpower)

        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the original Hodgkin Huxley equation.

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        # delta_Q = -(dyna.maxDmNa*P*(V - self.vrev))

        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem = sim.cc_cells[sim.iNa][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env = sim.cc_env[sim.iNa][cells.map_mem2ecm]

        else:
            c_env = sim.cc_env[sim.iNa]

        IdM = np.ones(sim.mdl)

        z_ion = sim.zs[sim.iNa] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmNa*P*1.0e-9

        # calculate specific ion flux contribution for this channel:
        delta_Q = stb.electroflux(c_env, c_mem, Dchan, p.tm * IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.clip_flux(delta_Q, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q, dyna.targets_vgNa, sim, cells, p)


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
class Nav1p2(VgNaABC):
    '''
    NaV model from Hammil et al 1991, derived from rat neocortical neurons.

    This channel produces well behaved action-potentials with a variety of vgK channels. Good general-purpose
    vgNa channel.

    Reference: Hammil, OP et al. Patch-clamp studies of voltage-gated currents in identified neurons
    of the rat cerebral cortex. Cereb. Cortex, 1991, 1, 48-61

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgNa channel type: Nav1.2')


        self.vrev = 50     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3  # FIXME 2 or 3 or user option?
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = (1 / (mAlpha + mBeta))/self.qt
        self._hInf = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))
        self._hTau = (1 / (
            (0.024 * ((V - 10.0) - -50.0)) /
            (1 - (np.exp(-((V - 10.0) - -50.0) / 5))) + (
                0.0091 * (-(V - 10.0) - 75.000123)) / (1 - (np.exp(-(-(V - 10) - 75.000123) / 5)))))/self.qt

class Nav1p3(VgNaABC):

    """
    Nav1.3 sodium channel TTX sensitive expressed in embryonic and neonatal tissue, rare in
    adults.

    It is thought that the fast activation and inactivation kinetics of Nav1.3, together with its rapid repriming
    kinetics and persistent current component, contributes to the development of spontaneous ectopic
    discharges and sustained rates of firing characteristics of injured sensory nerves

    Nav1.3 is predominantly expressed during pre-myelination stages of development in somatodendritic and axonal
    compartments of embryonic neurons. In rodents, the expression is attenuated after birth while in humans,
    the channel is highly expressed in somatodendritic compartment of myelinated neurons. Axotomy or other forms of
    nerve damage lead to the reexpression of NaV1.3 and the associated beta-3 subunit in sensory neurons,
    but not in primary motor neurons

    reference: Cummins et al. Nav1.3 sodium channel: Rapid repriming and slow closed-state
    inactivation display quantitative differences after expression in a mammalian
    cell line and in spinal sensory neurons. J. Neurosci., 2001, 21, 5952-61

    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel: Nav1p3')

        self.vrev = 50     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * ((V) - -26)) / (1 - (np.exp(-((V) - -26) / 9)))

        mBeta = (0.124 * (-(V) - 26)) / (1 - (np.exp(-(-(V) - 26) / 9)))

        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1 / (1 + np.exp((V - (-65.0)) / 8.1))

    def _calculate_state(self, V, dyna, sim, p):


        mAlpha = (0.182 * ((V) - -26)) / (1 - (np.exp(-((V) - -26) / 9)))

        mBeta = (0.124 * (-(V) - 26)) / (1 - (np.exp(-(-(V) - 26) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf = 1 / (1 + np.exp((V - (-65.0)) / 8.1))
        self._hTau = 0.40 + (0.265 * np.exp(-V / 9.47))

class NavRat2(VgNaABC):

    """
    Generic vgNa channel, slow to inactivate.

    Reference:  	McCormick DA. et al. A model of the electrophysiological properties of
    thalamocortical relay neurons. J. Neurophysiol., 1992 Oct , 68 (1384-400).
    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel: NavRat2')

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = 0.091 * (V + 38) / (1 - np.exp((-V - 38) / 5))
        mBeta = -0.062 * (V + 38) / (1 - np.exp((V + 38) / 5))

        dyna.m_Na = mAlpha / (mAlpha + mBeta)

        hAlpha = 0.016 * np.exp((-55 - V) / 15)
        hBeta = 2.07 / (np.exp((17 - V) / 21) + 1)

        dyna.h_Na = hAlpha / (hAlpha + hBeta)


    def _calculate_state(self, V, dyna, sim, p):

        mAlpha = 0.091 * (V + 38) / (1 - np.exp((-V - 38) / 5))
        mBeta = -0.062 * (V + 38) / (1 - np.exp((V + 38) / 5))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        hAlpha = 0.016 * np.exp((-55 - V) / 15)
        hBeta = 2.07 / (np.exp((17 - V) / 21) + 1)

        self._hInf = hAlpha / (hAlpha + hBeta)
        self._hTau = 1 / (hAlpha + hBeta)

class NavRat1(VgNaABC):

    """
    Generic vgNa channel, good activity.

    Reference: Huguenard JR. et al. Developmental changes in Na+ conductances in rat neocortical neurons:
     appearance of a slowly inactivating component. J. Neurophysiol., 1988 Mar , 59 (778-95).
    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel: NavRat1')

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1.0 / (1 + np.exp((V - -65) / 6.2))


    def _calculate_state(self, V, dyna, sim, p):


        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        self._hInf = 1.0 / (1 + np.exp((V - -65) / 6.2))
        self._hTau = 1 / ((0.024 * (V - -50)) / (1 - (np.exp(-(V - -50) / 5))) + (0.0091 * (-V - 75.000123)) / (
        1 - (np.exp(-(-V - 75.000123) / 5))))

class NavRat3(VgNaABC):  # FIXME finish this!

    """
    Generic vgNa channel.

    Reference:  	McCormick DA. et al. A model of the electrophysiological properties of
    thalamocortical relay neurons. J. Neurophysiol., 1992 Oct , 68 (1384-400).
    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel: NavRat3')

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1.0 / (1 + np.exp((V - -65) / 6.2))


    def _calculate_state(self, V, dyna, sim, p):


        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        self._hInf = 1.0 / (1 + np.exp((V - -65) / 6.2))
        self._hTau = 1 / ((0.024 * (V - -50)) / (1 - (np.exp(-(V - -50) / 5))) + (0.0091 * (-V - 75.000123)) / (
        1 - (np.exp(-(-V - 75.000123) / 5))))

class NaLeak(VgNaABC):

    '''
    Simple sodium leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a substance-modulated Na+ channel')


        self.vrev = 50     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        self.qt = 1.0

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Na = 1
        dyna.h_Na = 1

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








