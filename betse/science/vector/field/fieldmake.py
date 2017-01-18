#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector field factories, producing instances of the :class:`VectorFieldCells`
class.
'''

# ....................{ IMPORTS                            }....................
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

# ....................{ MAKERS ~ currents                  }....................
@type_check
def make_currents_intra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorFieldCells:
    '''
    Create and return a vector field cache of the current densities of all
    intracellular spaces for all time steps of the passed simulation, originally
    spatially situated at cell centres.

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
        Intracellular vector field cache of all current densities.
    '''

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
    intra- and extracellular spaces for all time steps of the passed simulation,
    originally spatially situated at environmental grid space centres.

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
        Intra- and extracellular vector field cache of all current densities.

    Raises
    ----------
    BetseSimConfigException
        If this simulation has disabled extracellular spaces.
    '''

    # If extracellular spaces are disabled, raise an exception.
    p.die_unless_ecm()

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(cells=cells, p=p, times_grids_centre=sim.I_tot_x_time),
        y=VectorCells(cells=cells, p=p, times_grids_centre=sim.I_tot_y_time),
        magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
    )

# ....................{ MAKERS ~ electric                  }....................
@type_check
def make_electric_extra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorFieldCells:
    '''
    Create and return a vector field cache of the eletric field across all
    extracellular spaces for all time steps of the passed simulation, originally
    spatially situated at environmental grid space centres.

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
        Extracellular electric vector field cache.

    Raises
    ----------
    BetseSimConfigException
        If this simulation has disabled extracellular spaces.
    '''

    # If extracellular spaces are disabled, raise an exception.
    p.die_unless_ecm()

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(
            cells=cells, p=p, times_grids_centre=sim.efield_ecm_x_time),
        y=VectorCells(
            cells=cells, p=p, times_grids_centre=sim.efield_ecm_y_time),
    )


@type_check
def make_electric_intra(
    sim:   'betse.science.sim.Simulator',
    cells: 'betse.science.cells.Cells',
    p:     'betse.science.parameters.Parameters',
) -> VectorFieldCells:
    '''
    Create and return a vector field cache of the eletric field across all
    intracellular spaces for all time steps of the passed simulation, originally
    spatially situated at cell membrane midpoints.

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
        Intracellular electric vector field cache.
    '''

    return VectorFieldCells(
        x=VectorCells(
            cells=cells, p=p, times_membranes_midpoint=sim.efield_gj_x_time),
        y=VectorCells(
            cells=cells, p=p, times_membranes_midpoint=sim.efield_gj_y_time),
    )
