#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Voltage-gated channel classes based on the Morris-Lecar simplified Hodgkin-Huxley channel approximations.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science.channels.channelsabc import ChannelsABC
from betse.util.io.log import logs

# .................... BASE                               ....................
class VgKABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all Voltage-gated Morris-Lecar channel classes.

    Attributes
    ----------
    _mInf :
        Equillibrium value function for m-gates of channel
    _mTau :
        Time-constant function for m-gates of channel

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

        if self.targets is None:

            V = vm * 1000

        else:
            V = vm[self.targets] * 1000

        self._init_state(V)

    def run(self, vm, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''

        if self.targets is None:

            V = vm*1000

        else:

            V = vm[self.targets] * 1000

        self._calculate_state(V)

        if self.kinetic_gate: # If the gate dynamics are time-dependent:

            self.update_ml(p, time_unit = self.time_unit)

        else: # the probability is simply the steady-state voltage-sensitivity:

            self.m = self._mInf

        if self.targets is None:

            self.P = self.m

        else:
            self.P = np.zeros(self.mdl)
            self.P[self.targets] = self.m


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

class Kv_ML1(VgKABC):
    '''
    Default delayed-rectifier voltage gated K+ channel model described in the Morris-Lecar literature.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kv_MorrisLecar1 ')

        self.time_unit = 1.0e3

        self.vrev = -85     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = True
        self.Phi = 0.066  # time-constant

        # initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
        self.m = 0.5*(1 + np.tanh((V - 12.0)/17.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V- 12)/17))
        self._mTau = 1/(np.cosh((V - 12)/(2*17)))

class Kv2p1_ML(VgKABC):
    '''
    Default delayed-rectifier voltage gated K+ channel model paralleling that of Korngreen et al, but using
    simplified Morris-Lecar formalism.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kv2.1, MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = True
        self.Phi = 0.066

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V - 14.0)/30.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V- 14)/30))
        self._mTau = 1/(np.cosh((V - 14)/(2*30)))

class Kv1p3_ML(VgKABC):
    '''
    Default delayed-rectifier voltage gated K+ channel model paralleling the Kv1.3 channel of  et al, but using
    simplified Morris-Lecar formalism.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kv1.3, Morris-Leclar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = True
        self.Phi = 0.066

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 14.0)/20.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 14)/20))
        self._mTau = 1/(np.cosh((V + 14)/(2*20)))

class Kv1p5_ML(VgKABC):
    '''
    Default delayed-rectifier voltage gated K+ channel model paralleling the Kv1.5 channel of  et al, but using
    simplified Morris-Lecar formalism.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kv_1.5, Morris-Lecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = True
        self.Phi = 0.066

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 6.0)/15.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 6)/15))
        self._mTau = 1/(np.cosh((V + 6)/(2*15)))

class Kv1p5S_ML(VgKABC):
    '''
    Slow version of delayed-rectifier voltage gated K+ channel model paralleling the Kv1.5 channel of  et al, but using
    simplified Morris-Lecar formalism, with extended time constant for large time-step simulations.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kv_1.5, Morris-Lecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = True
        self.Phi = 0.00066*2
        # self.Phi = 0.00066

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 6.0)/15.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 6)/15))
        self._mTau = 1/(np.cosh((V + 6)/(2*15)))

class Nav_ML(VgKABC):
    '''
    Persistent (i.e. inactivating) NaV channel current based on NaV1.6.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the NaV channel: Nav_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = 50     # reversal voltage used in model [mV]

        self.ions = ['Na']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 17.0)/18.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 17)/18))
        self._mTau = 1

class Cav_L_ML(VgKABC):
    '''
    Persistent (i.e. inactivating) L-type CaV channel current based on model of Carlin et al (a CaV1.3 channel).

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the CaV channel: Cav_L_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = 80     # reversal voltage used in model [mV]

        self.ions = ['Ca']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 20.0)/24.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 20)/24))
        self._mTau = 1

class Cav_L_ML2(VgKABC):
    '''
    Persistent (i.e. inactivating) L-type CaV channel current based on model of Carlin et al (a CaV1.3 channel).

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the CaV channel: Cav_L_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = 80     # reversal voltage used in model [mV]

        self.ions = ['Ca']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 20.0)/12.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 20)/12))
        self._mTau = 1

class Cav_N_ML(VgKABC):
    '''
    Persistent (i.e. inactivating) N-type CaV channel current based on CaV2.2 model of Huang et al.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the CaV channel: Cav_N_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = 80     # reversal voltage used in model [mV]

        self.ions = ['Ca']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 1.0)/18.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 1.0)/18.0))
        self._mTau = 1

class Cav_T_ML(VgKABC):
    '''
    Persistent (i.e. inactivating) T-type CaV channel current based on CaV3.1 model of Traboulsie et al.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the CaV channel: Cav_T_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = 80     # reversal voltage used in model [mV]

        self.ions = ['Ca']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 43.0)/24.0))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 43)/24.0))
        self._mTau = 1

class Kir_ML(VgKABC):
    '''
    Inward-rectifier Kir voltage-sensitive K+ channel current based on Kir2.1 model of Samy et al.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the KV channel: Kir_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K']
        self.rel_perm = [1]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 1/(np.cosh((V + 135)/(37)))

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 1/(np.cosh((V + 135)/(37)))
        self._mTau = 1


class HCN2_ML(VgKABC):
    '''
    HCN2 cation channel current based on HCN2 model of Moosmang et al.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the channel: HCN2_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K', 'Na']
        self.rel_perm = [1, 0.2]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 99.0)/12.4))  #19.2

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 99.0)/12.4))
        self._mTau = 1

class HCN4_ML(VgKABC):
    '''
    HCN4 cation channel current based on HCN2 model of Moosmang et al.

    '''
    def _init_state(self, V):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the channel: HCN2_MorrisLecar ')

        self.time_unit = 1.0e3

        self.vrev = -65     # reversal voltage used in model [mV]

        self.ions = ['K', 'Na']
        self.rel_perm = [1, 0.2]
        self.kinetic_gate = False
        self.Phi = 1.0

        # initialize values of the m gate of the channel based on m_inf:
        self.m = 0.5*(1 + np.tanh((V + 99.0)/19.2))  #19.2

    def _calculate_state(self, V):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self._mInf = 0.5*(1 + np.tanh((V + 99.0)/19.2))
        self._mTau = 1
