#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
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

    def init(self, dyna, sim, p):
        '''
        Initialize targeted voltage-gated sodium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''
        self.v_corr = 0  # correction factor for voltages

        V = sim.vm[dyna.targets_vgNa] * 1000 + self.v_corr

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
            V=sim.vm[dyna.targets_vgNa] * 1000 + self.v_corr,
            dyna=dyna, sim=sim, p=p)

        self._implement_state(dyna,sim,p)

    def _implement_state(self, dyna, sim, p):

        # calculate m and h channel states using RK4:
        dmNa = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhNa = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Na = dmNa(dyna.m_Na, p.dt*1e3) + dyna.m_Na
        dyna.h_Na = dhNa(dyna.h_Na, p.dt*1e3) + dyna.h_Na

        # calculate the open-probability of the channel:
        P = (dyna.m_Na ** self._mpower) * (dyna.h_Na ** self._hpower)

        # print(P.max())

        # # Find inds that should be zero, but are "hanging" the system:
        # Pcutoff = (P<0.0001).nonzero()
        # dyna.m_Na[Pcutoff] = self._mInf[Pcutoff]
        # dyna.h_Na[Pcutoff] = self._hInf[Pcutoff]

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa * P


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

    Note, the model has been adapted to have m-power =2 instead of 3, as 3 was found to not work well to produce
    action-potentials with the BETSE model.

    Reference: Hammil, OP et al. Patch-clamp studies of voltage-gated currents in identified neurons
    of the rat cerebral cortex. Cereb. Cortex, 1991, 1, 48-61

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgNa channel type: Nav1.2')

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
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))
        self._hTau = 1 / (
            (0.024 * ((V - 10.0) - -50.0)) /
            (1 - (np.exp(-((V - 10.0) - -50.0) / 5))) + (
                0.0091 * (-(V - 10.0) - 75.000123)) / (1 - (np.exp(-(-(V - 10) - 75.000123) / 5))))


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





