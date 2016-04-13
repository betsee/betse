#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated sodium channel classes.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.tissue.channels.channels_abc import ChannelsABC
from abc import ABCMeta, abstractmethod
from betse.util.io import logs

# ....................{ BASE                               }....................
class VgNaABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated sodium channel classes.

    Attributes
    ----------
    _mInf :
        Document me.
    _mTau :
        Document me.
    _hInf :
        Document me.
    _hTau :
        Document me.
    _v_inds_h :
        Document me.
    _v_inds_m :
        Document me.
    _corr_const :
        Document me.
    '''

    def init(self, dyna, sim, p):
        '''
        Initialize targeted voltage-gated sodium channels at the initial
        time step of the simulation based on the initial cell Vmems.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.
        '''

        V = sim.vm[dyna.targets_vgNa] * 1000

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
            V=sim.vm[dyna.targets_vgNa] * 1000,
            dyna=dyna, sim=sim, p=p)




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
class vgNa_Default(VgNaABC):
    '''
    NaV model from Hammil et al 1991, rat neocortical neurons.

    This channel produces well behaved action-potentials with a variety of vgK channels. Good general-purpose
    vgNa channel.

    Reference: Hammil, OP et al. Patch-clamp studies of voltage-gated currents in identified neurons
    of the rat cerebral cortex. Cereb. Cortex, 1991, 1, 48-61

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgNa channel type: vgNa_Default')
        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -39
        truth_inds_hb = V > -41

        v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -24
        truth_inds_mb = V > -26

        v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        # define small correction constant on the voltage
        self._corr_const = 1.0e-6

        # bump the V at nulclines with a small perturbation factor
        V[v_inds_m] = V + self._corr_const
        V[v_inds_h] = V + self._corr_const

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """
        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -39
        truth_inds_hb = V > -41

        self._v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -24
        truth_inds_mb = V > -26

        self._v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        # small correction constant on the voltage
        # self._corr_const = 1.0e-6

        V[self._v_inds_m] = V + self._corr_const
        V[self._v_inds_h] = V + self._corr_const

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))
        self._hTau = 1 / (
            (0.024 * ((V - 10.0) - -50.0)) /
            (1 - (np.exp(-((V - 10.0) - -50.0) / 5))) + (
                0.0091 * (-(V - 10.0) - 75.000123)) / (1 - (np.exp(-(-(V - 10) - 75.000123) / 5))))

        # calculate m and h channel states:
        dyna.m_Na = ((self._mInf - dyna.m_Na) / self._mTau) * p.dt * 1e3 + dyna.m_Na
        dyna.h_Na = ((self._hInf - dyna.h_Na) / self._hTau) * p.dt * 1e3 + dyna.h_Na

        # correction strategy for voltages (to prevent stalling at nullclines)----------
        dyna.m_Na[self._v_inds_m] = dyna.m_Na * (1 + self._corr_const)
        dyna.m_Na[self._v_inds_h] = dyna.m_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_m] = dyna.h_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_h] = dyna.h_Na * (1 + self._corr_const)

        # find indices where m or h are above or below the allowed 0 to 1 probability threshold:
        inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
        inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
        inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
        inds_hNa_over = (dyna.h_Na > 1.0).nonzero()

        dyna.m_Na[inds_mNa_over] = 1.0
        dyna.m_Na[inds_mNa_under] = 0.0
        dyna.h_Na[inds_hNa_over] = 1.0
        dyna.h_Na[inds_hNa_under] = 0.0

        # calculate the open-probability of the channel:
        P = (dyna.m_Na ** self._mpower) * (dyna.h_Na ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa * P


class vgNa_Rat(VgNaABC):

    def _init_state(self, V, dyna, sim, p):
        pass

    def _calculate_state(self, V, dyna, sim, p):
        pass


class vgNa_Squid(VgNaABC):
    '''
    Original Hodgkin-Huxley squid axon voltage gated sodium channel model.

    Reference: Hodgkin, AL et al. A quantitative description of membrane current and its application
    to conduction and excitation in nerve. 1952. Bull. Math. Biol., 1990, 52, 25-71.

    '''

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel type: vgNa_Squid')

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        # perturbation factor (prevents stalling)
        self._corr_const = 1.0e-6

        mAlpha = (0.1 * (25 - V)) / (np.exp((25 - V) / 10) - 1.0)
        mBeta = 4.0 * (np.exp(-V / 18))
        dyna.m_Na = mAlpha / (mAlpha + mBeta)

        hAlpha = 0.07 * np.exp(-V / 20)
        hBeta = 1 / (np.exp((30 - V) / 10) + 1.0)
        dyna.h_Na = hAlpha / (hAlpha + hBeta)


    def _calculate_state(self, V, dyna, sim, p):

        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -24
        truth_inds_hb = V > -26

        self._v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -24
        truth_inds_mb = V > -26
        #
        self._v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        V[self._v_inds_m] = V + self._corr_const
        V[self._v_inds_h] = V + self._corr_const

        mAlpha = (0.1 * (25 - V)) / (np.exp((25 - V) / 10) - 1.0)
        mBeta = 4.0 * (np.exp(-V / 18))
        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        hAlpha = 0.07 * np.exp(-V / 20)
        hBeta = 1 / (np.exp((30 - V) / 10) + 1.0)
        self._hInf = hAlpha / (hAlpha + hBeta)
        self._hTau = 1 / (hAlpha + hBeta)

        # calculate m and h channel states:
        dyna.m_Na = ((self._mInf - dyna.m_Na) / self._mTau) * p.dt * 1e3 + dyna.m_Na
        dyna.h_Na = ((self._hInf - dyna.h_Na) / self._hTau) * p.dt * 1e3 + dyna.h_Na

        # correction strategy for voltages (to prevent stalling at nullclines)----------
        dyna.m_Na[self._v_inds_m] = dyna.m_Na * (1 + self._corr_const)
        dyna.m_Na[self._v_inds_h] = dyna.m_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_m] = dyna.h_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_h] = dyna.h_Na * (1 + self._corr_const)

        # find indices where m or h are above or below the allowed 0 to 1 probability threshold:
        inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
        inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
        inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
        inds_hNa_over = (dyna.h_Na > 1.0).nonzero()

        dyna.m_Na[inds_mNa_over] = 1.0
        dyna.m_Na[inds_mNa_under] = 0.0
        dyna.h_Na[inds_hNa_over] = 1.0
        dyna.h_Na[inds_hNa_under] = 0.0

        # calculate the open-probability of the channel:
        P = (dyna.m_Na ** self._mpower) * (dyna.h_Na ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa * P

class vgNa_1p3(VgNaABC):

    """
    Nav1.3 sodium channel TTX sensitive expressed in embryonic and neonatal tissue, rare in
    adults.

    reference: Cummins et al. Nav1.3 sodium channel: Rapid repriming and slow closed-state
    inactivation display quantitative differences after expression in a mammalian
    cell line and in spinal sensory neurons. J. Neurosci., 2001, 21, 5952-61
    """

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel type: vgNa_1p3')

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        # perturbation factor (prevents stalling)
        self._corr_const = 1.0e-6

        mAlpha = (0.182 * ((V) - -26)) / (1 - (np.exp(-((V) - -26) / 9)))

        mBeta = (0.124 * (-(V) - 26)) / (1 - (np.exp(-(-(V) - 26) / 9)))

        dyna.m_Na = mAlpha / (mAlpha + mBeta)
        dyna.h_Na = 1 / (1 + np.exp((V - (-65.0)) / 8.1))

    def _calculate_state(self, V, dyna, sim, p):

        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -25
        truth_inds_hb = V > -27

        self._v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -25
        truth_inds_mb = V > -27
        #
        self._v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        V[self._v_inds_m] = V + self._corr_const
        V[self._v_inds_h] = V + self._corr_const


        mAlpha = (0.182 * ((V) - -26)) / (1 - (np.exp(-((V) - -26) / 9)))

        mBeta = (0.124 * (-(V) - 26)) / (1 - (np.exp(-(-(V) - 26) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)
        self._hInf = 1 / (1 + np.exp((V - (-65.0)) / 8.1))
        self._hTau = 0.40 + (0.265 * np.exp(-V / 9.47))

        # calculate m and h channel states:
        dyna.m_Na = ((self._mInf - dyna.m_Na) / self._mTau) * p.dt * 1e3 + dyna.m_Na
        dyna.h_Na = ((self._hInf - dyna.h_Na) / self._hTau) * p.dt * 1e3 + dyna.h_Na

        # correction strategy for voltages (to prevent stalling at nullclines)----------
        dyna.m_Na[self._v_inds_m] = dyna.m_Na * (1 + self._corr_const)
        dyna.m_Na[self._v_inds_h] = dyna.m_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_m] = dyna.h_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_h] = dyna.h_Na * (1 + self._corr_const)

        # find indices where m or h are above or below the allowed 0 to 1 probability threshold:
        inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
        inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
        inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
        inds_hNa_over = (dyna.h_Na > 1.0).nonzero()

        dyna.m_Na[inds_mNa_over] = 1.0
        dyna.m_Na[inds_mNa_under] = 0.0
        dyna.h_Na[inds_hNa_over] = 1.0
        dyna.h_Na[inds_hNa_under] = 0.0

        # calculate the open-probability of the channel:
        P = (dyna.m_Na ** self._mpower) * (dyna.h_Na ** self._hpower)

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa * P

class vgNa_1p6(VgNaABC):

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the vgNa channel type: vgNa_1p6')

        self._mpower = 1.0
        self._hpower = 0.0

        self._corr_const  = 1.0e-6

        dyna.m_Na = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        dyna.h_Na = 1.0

    def _calculate_state(self, V, dyna, sim, p):

        self._v_inds_m = []
        self._v_inds_h = []

        self._mInf = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self._mTau = 1

        self._hInf = 1.0
        self._hTau = 1.0

        # calculate m and h channel states:
        dyna.m_Na = ((self._mInf - dyna.m_Na) / self._mTau) * p.dt * 1e3 + dyna.m_Na
        dyna.h_Na = 1

        # correction strategy for voltages (to prevent stalling at nullclines)----------
        dyna.m_Na[self._v_inds_m] = dyna.m_Na * (1 + self._corr_const)
        dyna.m_Na[self._v_inds_h] = dyna.m_Na * (1 + self._corr_const)

        # find indices where m or h are above or below the allowed 0 to 1 probability threshold:
        inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
        inds_mNa_over = (dyna.m_Na > 1.0).nonzero()

        dyna.m_Na[inds_mNa_over] = 1.0
        dyna.m_Na[inds_mNa_under] = 0.0

        # calculate the open-probability of the channel:
        P = dyna.m_Na

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa * P



