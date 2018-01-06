#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: Consider moving this submodule into the channel-specific
#"betse.science.tissue.channels" subpackage, for disambiguity.

import numpy as np

# calcium (and voltage) gated K+ channel-------------------------------------------------------------------------------
def cagPotassium(dyna,sim,cells,p):
    """
    Model of a high-conductance calcium activated potassium channel (SK channel), obtained from
    the allosteric model of Cox DH, Cui J, Aldrich RW. J Gen Physiology. 1997. 110: 257-281.
    """

    # get data on cytosolic calcium levels in target cells:

    ca = sim.cc_cells[sim.iCa][cells.mem_to_cells][dyna.targets_cagK]

    # calculate different terms:
    t1 = 1 + ((ca*1e6)/10.22)
    t2 = 1 + ((ca*1e6)/0.89)

    # calculate probability of channel being open or closed:
    P = 1/(1+((t1/t2)**4)*6182*np.exp(-(1.64*p.F*sim.vm)/(p.R*sim.T)))

    # print(P.min(), P.mean(), P.max())

    # # ensure proper probability behaviour:
    # inds_P_over = (P > 1.0).nonzero()
    # P[inds_P_over] = 1.0
    #
    # inds_P_under = (P < 0.0).nonzero()
    # P[inds_P_under] = 0.0

    # calculate conductance of this potassium channel:
    # print(P.mean())

    # delta_Q = - (dyna.maxDmKcag * P * (sim.vm - self.vrev))
    sim.Dm_cag[sim.iK] = dyna.maxDmKcag*P
