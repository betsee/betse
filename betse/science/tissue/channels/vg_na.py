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

    def init(self, dyna, sim, cells, p):
        '''
        '''

        V = sim.vm[dyna.targets_vgNa] * 1000

        self._init_state(V=V, dyna=dyna, sim=sim, cells=cells, p=p)

        #FIXME: This looks totally bjorked.
        dyna.m_Na = self._mAlpha / (self._mAlpha + self._mBeta)
        dyna.h_Na = 1.0 / (1 + np.exp((V - -65.0 - 10.0) / 6.2))


    def run(self, dyna, sim, cells, p):
        '''
        Handle all **targeted voltage-gated sodium channels** (i.e., only
        applicable to specific tissue profiles) specified by the passed
        user-specified parameters on the passed tissue simulation and cellular
        world for the passed time step.

        Channel model uses Hodgkin-Huxley model for voltage gated sodium channels.
        '''

        self._calculate_state(
            V=sim.vm[dyna.targets_vgNa] * 1000,
            dyna=dyna, sim=sim, cells=cells, p=p)

        # calculate m and h channel states:
        dyna.m_Na = ((self._mInf - dyna.m_Na) / self._mTau) * p.dt * 1e3 + dyna.m_Na
        dyna.h_Na = ((self._hInf - dyna.h_Na) / self._hTau) * p.dt * 1e3 + dyna.h_Na

        # correction strategy for voltages (to prevent stalling at nullclines)----------
        dyna.m_Na[self._v_inds_m] = dyna.m_Na * (1 + self._corr_const)
        dyna.m_Na[self._v_inds_h] = dyna.m_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_m] = dyna.h_Na * (1 + self._corr_const)
        dyna.h_Na[self._v_inds_h] = dyna.h_Na * (1 + self._corr_const)


    @abstractmethod
    def _init_state(self, dyna, sim, cells, p):
        '''
        Do something.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, dyna, sim, cells, p):
        '''
        Do something.
        '''
        pass

# ....................{ SUBCLASS                           }....................
class VgNaHammil(VgNaABC):
    '''
    NaS model Hammil 1991

    This channel produces cool, well behaved AP with `.sim_ECM=True`. Slow to activate. Use!
    '''

    def _init_state(self, V, dyna, sim, cells, p):
        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -39
        truth_inds_hb = V > -41

        v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -24
        truth_inds_mb = V > -26

        v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        # small correction constant on the voltage
        corr_const = 1.0e-6

        V[v_inds_m] = V + corr_const
        V[v_inds_h] = V + corr_const

        self._mAlpha = (0.182 * ((V - 10.0) - -35.0)) / (1 - (np.exp(-((V - 10.0) - -35.0) / 9)))
        self._mBeta = (0.124 * (-(V - 10.0) - 35.0)) / (1 - (np.exp(-(-(V - 10.0) - 35.0) / 9)))


    def _calculate_state(self, V, dyna, sim, cells, p):
        # Find areas where the differential equation is intrinsically ill-behaved:
        truth_inds_ha = V < -39
        truth_inds_hb = V > -41

        self._v_inds_h = (truth_inds_ha * truth_inds_hb).nonzero()

        truth_inds_ma = V < -24
        truth_inds_mb = V > -26

        self._v_inds_m = (truth_inds_ma * truth_inds_mb).nonzero()

        # small correction constant on the voltage
        self._corr_const = 1.0e-6

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


class VgNaRat(VgNaABC):
    pass


class VgNaSquid(VgNaABC):
    '''
    HH squid model modified (sustains pulses).
    '''


    #FIXME: This method should be doing something.
    def _init_state(self, V, dyna, sim, cells, p):
        pass


    #FIXME: This method should also define:
    #"self._v_inds_m", "self._v_inds_h", "self._corr_const"
    def _calculate_state(self, V, dyna, sim, cells, p):
        V = V + 50

        mAlpha = (0.1*(25-V))/(np.exp((25-V)/10) -1.0)
        mBeta = 4.0 * (np.exp(-V/18))
        self._mInf = mAlpha/(mAlpha + mBeta)
        self._mTau = 1/(mAlpha + mBeta)
        hAlpha = 0.07 * np.exp(-V/20)
        hBeta = 1/(np.exp((30-V)/10) + 1.0)
        self._hInf = hAlpha/(hAlpha + hBeta)
        self._hTau = 1/(hAlpha + hBeta)

        # print(dyna.m_Na.min(),dyna.m_Na.max(),dyna.h_Na.min(),dyna.h_Na.max())