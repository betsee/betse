#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated sodium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channelo.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

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

    def init(self, vm, cells, p, targets = None):
        '''
        Initialize targeted voltage-gated sodium channels at the initial
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

        self.ions = ['Na']
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
class Nav1p2(VgNaABC):
    '''
    NaV model from Hammil et al 1991, derived from rat neocortical neurons.

    This channel produces well behaved action-potentials with a variety of vgK channels. Good general-purpose
    vgNa channel.

    Reference: Hammil, OP et al. Patch-clamp studies of voltage-gated currents in identified neurons
    of the rat cerebral cortex. Cereb. Cortex, 1991, 1, 48-61

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgNa channel type: Nav1.2')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0e3

        self.vrev = 50     # reversal voltage used in model [mV]

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = mAlpha / (mAlpha + mBeta)
        self.h = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = (1 / (mAlpha + mBeta))
        self._hInf = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))
        self._hTau = (1 / (
            (0.024 * ((V - 10.0) - -50.0)) /
            (1 - (np.exp(-((V - 10.0) - -50.0) / 5))) + (
                0.0091 * (-(V - 10.0) - 75.000123)) / (1 - (np.exp(-(-(V - 10) - 75.000123) / 5)))))

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

    def _init_state(self, V):

        logs.log_info('You are using the vgNa channel: Nav1p3')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0e3

        self.vrev = 50     # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * ((V) - -26)) / (1 - (np.exp(-((V) - -26) / 9)))

        mBeta = (0.124 * (-(V) - 26)) / (1 - (np.exp(-(-(V) - 26) / 9)))

        self.m = mAlpha / (mAlpha + mBeta)
        self.h = 1 / (1 + np.exp((V - (-65.0)) / 8.1))

    def _calculate_state(self, V):


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

    def _init_state(self, V):

        logs.log_info('You are using the vgNa channel: NavRat2')

        self.time_unit = 1.0e3

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = 0.091 * (V + 38) / (1 - np.exp((-V - 38) / 5))
        mBeta = -0.062 * (V + 38) / (1 - np.exp((V + 38) / 5))

        self.m = mAlpha / (mAlpha + mBeta)

        hAlpha = 0.016 * np.exp((-55 - V) / 15)
        hBeta = 2.07 / (np.exp((17 - V) / 21) + 1)

        self.h = hAlpha / (hAlpha + hBeta)


    def _calculate_state(self, V):

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

    def _init_state(self, V):

        logs.log_info('You are using the vgNa channel: NavRat1')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0e3

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        self.m = mAlpha / (mAlpha + mBeta)
        self.h = 1.0 / (1 + np.exp((V - -65) / 6.2))


    def _calculate_state(self, V):


        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        self._mInf = mAlpha / (mAlpha + mBeta)
        self._mTau = 1 / (mAlpha + mBeta)

        self._hInf = 1.0 / (1 + np.exp((V - -65) / 6.2))
        self._hTau = 1 / ((0.024 * (V - -50)) / (1 - (np.exp(-(V - -50) / 5))) + (0.0091 * (-V - 75.000123)) / (
        1 - (np.exp(-(-V - 75.000123) / 5))))

class NavRat3(VgNaABC):

    """
    Generic vgNa channel.

    Reference:  	McCormick DA. et al. A model of the electrophysiological properties of
    thalamocortical relay neurons. J. Neurophysiol., 1992 Oct , 68 (1384-400).
    """

    def _init_state(self, V):

        logs.log_info('You are using the vgNa channel: NavRat3')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0e3

        self.vrev = 50  # reversal voltage used in model [mV]

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 3
        self._hpower = 1

        mAlpha = (0.182 * (V - -35)) / (1 - (np.exp(-(V - -35) / 9)))
        mBeta = (0.124 * (-V - 35)) / (1 - (np.exp(-(-V - 35) / 9)))

        self.m = mAlpha / (mAlpha + mBeta)
        self.h = 1.0 / (1 + np.exp((V - -65) / 6.2))


    def _calculate_state(self, V):


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

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using an Na+ leak channel')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0

        self.vrev = 50     # reversal voltage used in model [mV]

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

class Nav1p6(VgNaABC):

    """
    Nav1.6 was detected during the embryonic period in brain. This is the most abundantly expressed
    NaV channel in the CNS during adulthood. Is the most abundant channel at mature nodes of Ranvier
    in myelinated axons in the CNS. NaV1.6 is broadly expressed in the nervous system in a variety
    of cells including Purkinje cells, motor neurons, pyramidal and granule neurons, glial cells and
    Schwann cells and is enriched at the nodes of Ranvier.
    Nav1.6 channels have been also detected in immune cells, such microglia and macrophagues
     and in cultured microglia, Nav1.6 is the most prominently expressed sodium channel.

     Reference: Smith MR. et al. Functional analysis of the mouse Scn8a sodium channel. \
     J. Neurosci., 1998 Aug 15 , 18 (6093-102).

    """

    def _init_state(self, V):

        logs.log_info('You are using the persistent vgNa channel: Nav1p6')

        # time scale of the model (most models are in miliseconds, some are in seconds)
        self.time_unit = 1.0e3

        self.vrev = 50     # reversal voltage used in model [mV]

        self._mpower = 1.0
        self._hpower = 0.0

        self.m = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self.h = 1.0

    def _calculate_state(self, V):

        self._mInf = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self._mTau = 1

        self._hInf = 1.0
        self._hTau = 1.0

# class Nav(VgNaABC):
#
#     """
#     Minimal, general model of voltage gated sodium channels, from Pospischil et al 2008.
#
#
#      Reference: Pospischil M. et al. Minimal Hodgkin–Huxley type models for different classes
#      of cortical and thalamic neurons. Biological Cybernetics., 2008.
#
#      This model was not obtained from Channelpedia.
#
#     """
#
#     def _init_state(self, V, sim, p):
#
#         logs.log_info('You are using the general Nav channel: Nav')
#
#         self.time_unit = 1.0
#
#         self.vrev = 50     # reversal voltage used in model [mV]
#
#         self._mpower = 3.0
#         self._hpower = 1.0
#
#
#         self._alpha_m = -(0.32*(V - 13))/(np.exp(-(V-13)/4)-1)
#         self._beta_m = -(0.28 * (V - 40)) / (np.exp((V - 40) / 5) - 1)
#
#         self._alpha_h = 0.128*np.exp((V - 17)/18)
#         self._beta_h = 4/ (np.exp(-(V - 40) / 5) + 1)
#
#         self.m = self._alpha_m/(self._alpha_m + self._beta_m)
#         self.h = self._alpha_h/(self._alpha_h + self._beta_h)
#
#     def _calculate_state(self, V, sim, p):
#
#         self._alpha_m = -(0.32*(V - 13))/(np.exp(-(V-13)/4)-1)
#         self._beta_m = -(0.28 * (V - 40)) / (np.exp((V - 40) / 5) - 1)
#
#         self._alpha_h = 0.128*np.exp((V - 17)/18)
#         self._beta_h = 4/ (np.exp(-(V - 40) / 5) + 1)











