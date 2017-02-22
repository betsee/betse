#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector factories, producing instances of the :class:`VectorCells` class.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.vector.vectorcls import VectorCells
from betse.util.type.types import type_check


# ....................{ MAKERS                             }....................
@type_check
def make_voltages_intra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorCells:
    '''
    Create and return a vector cache of the voltages across all intracellular
    membranes (i.e., gap junctions) for all time steps of the passed simulation.

    Parameters
    ----------
    sim : Simulator
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.

    Returns
    ----------
    VectorFieldCells
        Intracellular vector cache of all voltages.
    '''

    return VectorCells(
        cells=cells, p=p,
        times_membranes_midpoint=expmath.upscale_cell_data(sim.vm_time),
    )
