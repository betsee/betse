#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector field factories, producing instances of the :class:`VectorFieldCells`
class.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimConfigException
from betse.science.vector.vectorcls import VectorCells
from betse.science.vector.field.fieldcls import VectorFieldCells
from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
_MAGNITUDE_FACTOR_CURRENTS = 100
'''
Factor by which to multiply each magnitude of each vector in each vector field
of intra- and/or extracellular current densities, producing magnitude in units
of uA/cm^2.
'''

# ....................{ MAKERS                             }....................
@type_check
def make_currents_intra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorFieldCells:
    '''
    Create and return a vector field cache of the current densities of all
    intracellular spaces for all time steps of the passed simulation.

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
        Vector field cache of all intracellular spaces.
    '''

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(cells=cells, p=p, times_cells_centre=sim.I_cell_x_time),
        y=VectorCells(cells=cells, p=p, times_cells_centre=sim.I_cell_y_time),
        magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
    )


@type_check
def make_currents_intra_extra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorFieldCells:
    '''
    Create and return a vector field cache of the current densities of all
    intra- and extracellular spaces for all time steps of the passed simulation.

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
        Vector field cache of all intra- and extracellular spaces.

    Raises
    ----------
    BetseSimConfigException
        If this simulation disabled extracellular spaces.
    '''

    # If this simulation disabled extracellular spaces, raise an exception.
    if not p.sim_ECM:
        raise BetseSimConfigException(
            'Extracellular spaces disabled, but required by the requested '
            'total current density vector field.')

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(cells=cells, p=p, times_grids_centre=sim.I_tot_x_time),
        y=VectorCells(cells=cells, p=p, times_grids_centre=sim.I_tot_y_time),
        magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
    )
