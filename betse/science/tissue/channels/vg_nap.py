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
from betse.science import sim_toolbox as stb


# .................... BASE                               ....................
class VgNaPABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all persistent Voltage-gated sodium channel classes.

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

        self.v_corr = 10  # offset of voltages in the model -- experimental junction voltage [mV]

        V = sim.vm[dyna.targets_vgNaP] * 1000 + self.v_corr

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)


    def run(self, dyna, sim, cells, p):
        '''
        Handle all targeted voltage-gated sodium channels by working with the passed
        user-specified parameters on the tissue simulation and cellular
        world for a time step.

        Channel model uses Hodgkin-Huxley kinetic model style
        for voltage gated channels.

        '''
        V = sim.vm[dyna.targets_vgNaP] * 1000 + self.v_corr

        self._calculate_state(V, dyna, sim, p)

        self._implement_state(V, dyna,sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):

        # calculate m and h channel states using RK4:
        dmNaP = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhNaP = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_NaP = dmNaP(dyna.m_NaP, p.dt*1e3) + dyna.m_NaP
        dyna.h_NaP = dhNaP(dyna.h_NaP, p.dt*1e3) + dyna.h_NaP

        # calculate the open-probability of the channel:
        P = (dyna.m_NaP ** self._mpower) * (dyna.h_NaP ** self._hpower)

        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the original Hodgkin Huxley equation.

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):
        delta_Q = -(dyna.maxDmNaP * P * (V - self.vrev)) / cells.mem_sa[dyna.targets_vgNaP]

        # the cube power in the vgNa expression is rather difficult mathematically, but necessary
        # clip the unreasonably high portions of the Na+ flux, so as not to overload the system:
        inds_over = (delta_Q > 1.0e-4).nonzero()
        delta_Q[inds_over] = 1.0e-4

        # update the fluxes across the membrane to account for charge transfer from HH flux:
        sim.fluxes_mem[sim.iNa][dyna.targets_vgNaP] = delta_Q

        # update the concentrations of Na in cells and environment using HH flux delta_Q:

        # first in cells:
        sim.cc_mems[sim.iNa][dyna.targets_vgNaP] = \
            sim.cc_mems[sim.iNa][dyna.targets_vgNaP] + \
            delta_Q * (cells.mem_sa[dyna.targets_vgNaP] / cells.mem_vol[dyna.targets_vgNaP]) * p.dt

        if p.sim_ECM is False:

            # transfer charge directly to the environment:

            sim.cc_env[sim.iNa][dyna.targets_vgNaP] = \
                sim.cc_env[sim.iNa][dyna.targets_vgNaP] - \
                delta_Q * (cells.mem_sa[dyna.targets_vgNaP] / cells.mem_vol[dyna.targets_vgNaP]) * p.dt

            # assume auto-mixing of environmental concs
            sim.cc_env[sim.iNa][:] = sim.cc_env[sim.iNa].mean()

        else:

            flux_env = np.zeros(sim.edl)
            flux_env[cells.map_mem2ecm][dyna.targets_vgNaP] = -delta_Q

            # save values at the cluster boundary:
            bound_vals = flux_env[cells.ecm_bound_k]

            # set the values of the global environment to zero:
            flux_env[cells.inds_env] = 0

            # finally, ensure that the boundary values are restored:
            flux_env[cells.ecm_bound_k] = bound_vals

            # Now that we have a nice, neat interpolation of flux from cell membranes, multiply by the,
            # true membrane surface area in the square, and divide by the true ecm volume of the env grid square,
            # to get the mol/s change in concentration (divergence):
            delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol

            # update the concentrations:
            sim.cc_env[sim.iNa][:] = sim.cc_env[sim.iNa][:] + delta_env * p.dt

        # update the concentration intra-cellularly:
        sim.cc_mems[sim.iNa], sim.cc_cells[sim.iNa], _ = stb.update_intra(sim, cells, sim.cc_mems[sim.iNa],
            sim.cc_cells[sim.iNa], sim.D_free[sim.iNa], sim.zs[sim.iNa], p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        sim.update_V(cells, p)

        # Define ultimate activity of the vgNa channel:
        # sim.Dm_vg2[sim.iNa][dyna.targets_vgNaP] = dyna.maxDmNaP * P


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

class Nav1p6(VgNaPABC):

    def _init_state(self, V, dyna, sim, p):

        logs.log_info('You are using the persistent vgNa channel: Nav1p6')

        self.vrev = 50     # reversal voltage used in model [mV]
        Texpt = 23    # temperature of the model in degrees C
        simT = sim.T - 273   # model temperature in degrees C
        # self.qt = 2.3**((simT-Texpt)/10)
        self.qt = 1.0  # FIXME implement this!

        self._mpower = 1.0
        self._hpower = 0.0

        dyna.m_NaP = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        dyna.h_NaP = 1.0

    def _calculate_state(self, V, dyna, sim, p):

        self._mInf = 1.0000 / (1 + np.exp(-0.03937 * 4.2 * (V - -17.000)))
        self._mTau = 1

        self._hInf = 1.0
        self._hTau = 1.0