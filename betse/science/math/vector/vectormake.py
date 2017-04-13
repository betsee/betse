#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector factories, producing instances of the :class:`VectorCellsCache` class.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.simulate.simphase import SimPhase
from betse.util.type.types import type_check

# ....................{ MAKERS                             }....................
@type_check
def make_voltages_membrane(phase: SimPhase) -> VectorCellsCache:
    '''
    Vector caching all **transmembrane voltages** (i.e., voltages across all
    gap junctions connecting intracellular membranes) for all time steps of the
    passed simulation phase.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCellsCache
        Vector caching all transmembrane voltages.
    '''

    return VectorCellsCache(
        phase=phase,
        times_membranes_midpoint=expmath.upscale_cell_data(phase.sim.vm_time),
    )