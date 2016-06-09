#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.exceptions import BetseExceptionLambda
from betse.util.type import types
from betse.science import toolbox as tb

# calcium (and voltage) gated K+ channel-------------------------------------------------------------------------------
#
def cagPotassium(dyna,sim,cells,p):
    """
    Model of a high-conductance calcium activated potassium channel (SK channel), obtained from
    the allosteric model of Cox DH, Cui J, Aldrich RW. J Gen Physiology. 1997. 110: 257-281.

    """

    # get data on cytosolic calcium levels in target cells
    if p.sim_ECM is False:

        ca = sim.cc_mems[sim.iCa][dyna.targets_cagK]

    else:

        ca = sim.cc_mems[sim.iCa][cells.mem_to_cells][dyna.targets_cagK]

    # calculate different terms:
    t1 = 1 + ((ca*1e3)/10.22)
    t2 = 1 + ((ca*1e3)/0.89)

    # calculate probability of channel being open or closed:
    P = 1/(1+((t1/t2)**4)*6182*np.exp(-(1.64*p.F*sim.vm)/(p.R*sim.T)))

    # ensure proper probability behaviour:
    inds_P_over = (P > 1.0).nonzero()
    P[inds_P_over] = 1.0

    inds_P_under = (P < 0.0).nonzero()
    P[inds_P_under] = 0.0

    # calculate conductance of this potassium channel:
    sim.Dm_cag[sim.iK] = dyna.maxDmKcag*P







