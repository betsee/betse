#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all channel classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
import numpy as np
from betse.science import sim_toolbox as stb

# ....................{ BASE                               }....................
class ChannelsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all channel classes.

    Attributes
    ----------
    '''

    @abstractmethod
    def init(self, dyna, sim, cells, p):
        '''
        Runs the initialization sequence for a voltage gated ion channel.
        '''
        pass

    @abstractmethod
    def run(self, dyna, sim, cells, p):
        '''
        Runs the voltage gated ion channel.
        '''
        pass

    def update_charge(self, ion_index, delta_Qo, targets, sim, cells, p):

        """
        A general helper function to update charge in the cell and environment
        given a flux derived from the GHK flux equation.

        Parameters
        ----------------
        ion_index:  index of an ion in the sim module (i.e. sim.iNa, sim.iK, sim.iCl, etc)
        delta_Qo:    GHK flux component for the channel state
        targets:    Indices to the cell membrane targets for the channel (e.g. dyna.targets_vgNa)
        sim:        Instance of sim object
        cells:      Instance of cells object
        p:          Instances of params object

        """

        # multiply by the modulator:

        cells.cell_i = np.asarray(cells.cell_i)

        targets = np.asarray(targets) # convert targets into an array so we can index with it

        # derive the targeted expression for delta_Q:
        delta_Q = np.zeros(sim.mdl)
        delta_Q[targets] = delta_Qo[targets]


        # update charge in the cell and environment, assuming a trans-membrane flux occurs due to open channel state,
        # which is described by the GHK flux equation.

        # update the fluxes across the membrane to account for charge transfer from channel flux:
        sim.fluxes_mem[ion_index] = sim.fluxes_mem[ion_index] + delta_Q

        # update the concentrations of ion in cells and environment using GHK derived flux delta_Q:
        master_inds = cells.cell_i[cells.mem_to_cells][targets]

        ccell =  sim.cc_cells[ion_index]

        # take the divergence of the flux:
        divQ = np.dot(cells.M_sum_mems, delta_Q*cells.mem_sa)/cells.cell_vol

        # first in cells:
        sim.cc_cells[ion_index] = (ccell + divQ*p.dt)

        # next on membranes:
        # sim.cc_at_mem[ion_index] = (sim.cc_at_mem[ion_index] +
        #                                      delta_Q*(cells.mem_sa/cells.mem_vol)*p.dt)

        if p.is_ecm is False:

            # transfer charge directly to the environment:
            sim.cc_env[ion_index] = (
                sim.cc_env[ion_index] -
                delta_Q * (cells.mem_sa / cells.mem_vol) * p.dt)

            # assume auto-mixing of environmental concs
            sim.cc_env[ion_index][:] = sim.cc_env[ion_index].mean()

        else:

            # delta_env = np.dot(cells.M_divmap_mem2ecm, -delta_Q)
            delta_env = stb.div_env(-delta_Q, cells,p)

            # update the concentrations:
            sim.cc_env[ion_index]= sim.cc_env[ion_index] + delta_env * p.dt

        # recalculate the net, unbalanced charge and voltage in each cell:
        # sim.update_V(cells, p)

    def clip_flux(self, delta_Q, threshold = 1.0e-3):
        """
        Clips flux so that it remains within stable limits of the BETSE model
        for a reasonable time step.

        delta_Q:  Flux
        threshold: Flux is clipped to within +/- threshold

        """

        inds_over = (delta_Q > threshold).nonzero()
        delta_Q[inds_over] = threshold

        inds_under = (delta_Q < -threshold).nonzero()
        delta_Q[inds_under] = -threshold

        return delta_Q

