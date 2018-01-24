#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated potassium channel classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb
from betse.science.channelo.channelsabc import ChannelsABC
from betse.science.math import toolbox as tb
from betse.util.io.log import logs

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

        self.ions = ['K']
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

class Kv1p1(VgKABC):
    '''
    Kv1.1 model from Christie et al, cloned from rat brain, studied in xenopus.

    Kv1.1 are low-voltage activated (LVA) channels, expressed primarily in the central nervous system, and
    which open with small depolarizations at or below resting potential

    Reference: Christie MJ. et al. Expression of a cloned rat brain potassium channel in Xenopus oocytes.
    Science, 1989 Apr 14 , 244 (221-4).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1p1 ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -30.5000) / -11.3943))
        self.h = 1.0000 / (1 + np.exp((V - -30.0000) / 27.3943))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 2


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -30.5000) / -11.3943))
        self._mTau = 30.0000 / (1 + np.exp((V - -76.5600) / 26.1479))
        self._hInf = 1.0000 / (1 + np.exp((V - -30.0000) / 27.3943))
        self._hTau = 15000.0000 / (1 + np.exp((V - -160.5600) / -100.0000))

class Kv1p2(VgKABC):
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

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1.2')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
        self.h = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
        self._mTau = 150.0000/(1+ np.exp((V + 67.5600)/34.1479))
        self._hInf = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))
        self._hTau = 15000.0000/(1+ np.exp(-(V + 46.5600)/44.1479))

class Kv1p3(VgKABC):
    '''
    Kv1.3 model from Douglass et al 1990.

    Potassium voltage-gated channel Kv1.3 is a member of the shaker-related subfamily and belongs to the
    delayed rectifier class of channels, which allow nerve cells to efficiently re-polarize following an
    action potential. This channel is implicated in a host of non-neural activity, including autoimmune
    and inflammatory disorders, cell proliferation, and cancer development.

    reference:	Douglass J. et al. Characterization and functional expression of a rat genomic DNA
    clone encoding a lymphocyte potassium channel. J. Immunol., 1990 Jun 15 , 144 (4841-50).


    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1.3')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -14.1000) / -10.3000))
        self.h = 1.0000 / (1 + np.exp((V - -33.0000) / 3.7000))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -14.1000) / -10.3000))
        self._mTau = (-0.2840 * V) + 19.1600
        self._hInf = 1.0000 / (1 + np.exp((V - -33.0000) / 3.7000))
        self._hTau = (-13.7600 * V) + 1162.4000

class Kv1p4(VgKABC):
    '''
    Kv1.4 model from Stuhmer et al 1989.



    reference: W Stühmer et. al; EMBO J. 1989 Nov


    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1.4')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V + 21.7000) / -16.9000))
        self.h = 1.0000 / (1 + np.exp((V + 73.6000) / 12.8000))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V + 21.7000) / -16.9000))
        self._mTau = 3.0
        self._hInf = 1.0000 / (1 + np.exp((V + 73.6000) / 12.8000))
        self._hTau = 119.0

class Kv1p5(VgKABC):
    '''
    Kv1.5 model from

    This channel produces well behaved action-potentials with a variety of vgNa channels. Good general-purpose
    vgK channel.

    Kv1.5 underlies the cardiac ultra-rapid delayed rectifier potassium current (IKur) in humans.
    In addition, Kv1.5 channels are also expressed in many other organs, such as pulmonary arteries,
    brain, skeletal muscle, and have crucial function in cell cycle regulation.

    Reference: Philipson LH. et al. SequenceTypes and functional expression in Xenopus oocytes of a
    human insulinoma and islet potassium channel. Proc. Natl. Acad. Sci. U.S.A., 1991 Jan 1 , 88 (53-7).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1p5 ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1.0000 / (1 + np.exp((V - -6.0000) / -6.4000))
        self.h = 1.0000 / (1 + np.exp((V - -25.3000) / 3.5000))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1.0000 / (1 + np.exp((V - -6.0000) / -6.4000))
        self._mTau = (-0.1163 * V) + 8.3300
        self._hInf = 1.0000 / (1 + np.exp((V - -25.3000) / 3.5000))
        self._hTau = (-15.5000 * V) + 1620.0000

class Kv1p6(VgKABC):
    '''
    Kv1.6 model from Grupe et al. 1990.



    Reference: A Grupe et. al; EMBO J. 1990 Jun

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: Kv1p6 ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(((V - (-20.800)) / (-8.100))))
        self.h =  1 / (1 + np.exp(((V - (-22.000)) / (11.390))))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp(((V - (-20.800)) / (-8.100))))
        self._mTau = 30.000 / (1 + np.exp(((V - (-46.560)) / (44.140))))
        self._hInf = 1 / (1 + np.exp(((V - (-22.000)) / (11.390))))
        self._hTau = 5000.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))

class Kv2p1(VgKABC):
    """
    Delayed rectifier potassium channel found widespread through many tissue types.

    Reference: 	VanDongen AM. et al. Alteration and restoration of K+ channel function by deletions
    at the N- and C-termini. Neuron, 1990 Oct , 5 (433-43).

    """


    def _init_state(self, V):
        """
        Run initialization calculation for m and h gates of the channel at starting Vmem value.
        """
        logs.log_info('You are using the vgK channel: Kv2p1 ')

        self.time_unit = 1.0e3

        self.vrev = -65  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(((V - (-9.200)) / (-6.600))))
        self.h = 1 / (1 + np.exp(((V - (-19.000)) / (5.000))))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V):

        self._mInf = 1 / (1 + np.exp(((V - (-9.200)) / (-6.600))))
        self._mTau = 100.000 / (1 + np.exp(((V - (-46.560)) / (44.140))))
        self._hInf = 1 / (1 + np.exp(((V - (-19.000)) / (5.000))))
        self._hTau = 10000.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))

class Kv2p2(VgKABC):
    """
    Reference: F Schmalz et. al; Am. J. Physiol. 1998 May


    """


    def _init_state(self, V):
        """
        Run initialization calculation for m and h gates of the channel at starting Vmem value.
        """
        logs.log_info('You are using the vgK channel: Kv2p2 ')

        self.time_unit = 1.0e3

        self.vrev = -65  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(((V - (5.000)) / (-12.000))))
        self.h = 1 / (1 + np.exp(((V - (-16.300)) / (4.800))))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V):

        self._mInf = 1 / (1 + np.exp(((V - (5.000)) / (-12.000))))
        self._mTau = 130.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))
        self._hInf = 1 / (1 + np.exp(((V - (-16.300)) / (4.800))))
        self._hTau = 10000.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))

class Kv3p1(VgKABC):
    """
 	J Rettig et. al; EMBO J. 1992 Jul

    """


    def _init_state(self, V):
        """
        Run initialization calculation for m and h gates of the channel at starting Vmem value.
        """
        logs.log_info('You are using the vgK channel: Kv3p1 ')

        self.time_unit = 1.0e3

        self.vrev = -65  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(((V - (18.700)) / (-9.700))))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 0

    def _calculate_state(self, V):

        self._mInf = 1 / (1 + np.exp(((V - (18.700)) / (-9.700))))
        self._mTau = 20.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))
        self._hInf = 1
        self._hTau = 1

class Kv3p2(VgKABC):
    """
    R Hernandez-Pineda et. al; J. Neurophysiol. 1999 Sep

    """


    def _init_state(self, V):
        """
        Run initialization calculation for m and h gates of the channel at starting Vmem value.
        """
        logs.log_info('You are using the vgK channel: Kv3p2 ')

        self.time_unit = 1.0e3

        self.vrev = -65  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp((V - -0.373267) / -8.568187))
        self.h = 1

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 0

    def _calculate_state(self, V):

        self._mInf =  1 / (1 + np.exp((V - -0.373267) / -8.568187))
        self._mTau = 3.241643 + (19.106496 / (1 + np.exp((V - 19.220623) / 4.451533)))
        self._hInf = 1
        self._hTau = 1

class Kv3p3(VgKABC):


    """

     Kv3p3 is a rapidly inactivating, A-type inward/outward potassium current, implicated in rapid firing of
     neurons and is primarily found in Purkunjie neurons.

     Reference: Desai R. et al. Protein kinase C modulates inactivation of Kv3.3 channels.
     J. Biol. Chem., 2008 Aug 8 , 283 (22283-94).

    """

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the A-current K+ channel: Kv3p3')

        self.time_unit = 1.0e3

        # FIXME: IS THIS RIGHT?????
        self.vrev = 82.0  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp((V - 35) / -7.3))
        self.h =  0.25 + (0.75 / (1 + np.exp((V - (-28.293856)) / 29.385636)))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 2
        self._hpower = 1


    def _calculate_state(self, V):

        self._mInf = 1 / (1 + np.exp((V - 35) / -7.3))
        self._mTau = 0.676808 + (27.913114 / (1 + np.exp((V - 22.414149) / 9.704638)))
        self._hInf = 0.25 + (0.75 / (1 + np.exp((V - (-28.293856)) / 29.385636)))
        self._hTau = 199.786728 + (2776.119438 * np.exp(-V / 7.309565))

class Kv3p4(VgKABC):

    """

     Kv3p4 is a rapidly inactivating, A-type inward/outward potassium current.

     Reference: Vega-Saenz de Miera E. et al. Cloning of ShIII (Shaw-like) cDNAs encoding a novel
     high-voltage-activating, TEA-sensitive, type-A K+ channel. Proc. Biol. Sci., 1992 Apr 22 , 248 (9-18).

    """

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the K+ channel: Kv3p4')

        self.time_unit = 1.0e3

        self.vrev = -65.0  # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(((V - (-3.400)) / (-8.400))))
        self.h =  1 / (1 + np.exp(((V - (-53.320)) / (7.400))))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp(((V - (-3.400)) / (-8.400))))
        self._mTau = 10.000 / (1 + np.exp(((V - (4.440)) / (38.140))))
        self._hInf = 1 / (1 + np.exp(((V - (-53.320)) / (7.400))))
        self._hTau = 20000.000 / (1 + np.exp(((V - (-46.560)) / (-44.140))))

class K_Fast(VgKABC):
    '''
    "K Fast" model from Korngreen et al.

    Reference: Korngreen A. et al. Voltage-gated K+ channels in layer 5 neocortical pyramidal
    neurones from young rats: subtypes and gradients. J. Physiol. (Lond.), 2000 Jun 15 , 525 Pt 3 (621-39).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the vgK channel: K_Fast ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the potassium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp(-(V + 47) / 29))
        self.h = 1 / (1 + np.exp(-(V + 56) / -10))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """
        self._mInf = 1 / (1 + np.exp(-(V + 47) / 29))
        self._mTau = (0.34 + 0.92 * np.exp(-((V + 71) / 59)**2))
        self._hInf = 1 / (1 + np.exp(-(V + 56) / -10))
        self._hTau = (8 + 49 * np.exp(-((V + 73) / 23)**2))

class KLeak(VgKABC):

    '''
    Simple potassium leak channel -- always open -- for substance modulation.

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using a K+ Leak channel')

        self.time_unit = 1.0

        self.vrev = -65     # reversal voltage used in model [mV]

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



        # print(self.modulator.mean())

class Kir2p1(VgKABC):
    '''
    Kir 2.1 model from Makary et al.

    This channel has a greater tendency to allow potassium to flow into a cell rather than out of a cell,
    probably participates in establishing action potential waveform and excitability of neuronal and muscle tissues.
    Mutations in this gene have been associated with Andersen syndrome, which is characterized by periodic paralysis,
    cardiac arrhythmias, and dysmorphic features

    Reference: Makary SM. et al. A difference in inward rectification and polyamine block and permeation between
     the Kir2.1 and Kir3.1/Kir3.4 K+ channels. J. Physiol. (Lond.), 2005 Nov 1 , 568 (749-66).

    '''

    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the inward rectifying K+ channel: Kir2p1')

        self.time_unit = 1.0e3

        self.vrev = -70.6     # reversal voltage used in model [mV]

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 1 / (1 + np.exp((V - (-96.48)) / 23.26))
        self.h = 1 / (1 + np.exp((V - (-168.28)) / -44.13))

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 1
        self._hpower = 2

        # self.rectification = -1


    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1 / (1 + np.exp((V - (-96.48)) / 23.26))
        self._mTau = 3.7 + (-3.37 / (1 + np.exp((V - -32.9) / 27.93)))
        self._hInf = 1 / (1 + np.exp((V - (-168.28)) / -44.13))
        self._hTau = 0.85 + (306.3 / (1 + np.exp((V - -118.29) / -27.23)))












