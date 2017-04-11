#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector field factories, producing instances of the :class:`VectorFieldCells`
class.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.vectorcls import VectorCells
from betse.science.math.vector.field.fieldcls import VectorFieldCells
from betse.science.simulate.simphase import SimPhase
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
def make_currents_intra(phase: SimPhase) -> VectorFieldCells:
    '''
    Vector field caching all intracellular current densities for all time steps
    of the passed simulation phase, originally spatially situated at cell
    centres.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCells
        Vector field caching all intracellular current densities.
    '''

    return VectorFieldCells(
        x=VectorCells(phase=phase, times_cells_centre=phase.sim.I_cell_x_time),
        y=VectorCells(phase=phase, times_cells_centre=phase.sim.I_cell_y_time),
        magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
    )


@type_check
def make_currents_extra(phase: SimPhase) -> VectorFieldCells:
    '''
    Vector field caching all exracellular current densities for all time steps
    of the passed simulation phase, originally spatially situated at
    environmental grid space centres.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCells
        Vector field caching all extracellular current densities.

    Raises
    ----------
    BetseSimConfigException
        If this simulation has disabled extracellular spaces.
    '''

    # If extracellular spaces are disabled, raise an exception.
    phase.p.die_unless_ecm()

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(phase=phase, times_grids_centre=phase.sim.I_tot_x_time),
        y=VectorCells(phase=phase, times_grids_centre=phase.sim.I_tot_y_time),
        magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
    )

# ....................{ MAKERS ~ electric                  }....................
@type_check
def make_electric_intra(phase: SimPhase) -> VectorFieldCells:
    '''
    Vector field caching the electric field across all intracellular spaces for
    all time steps of the passed simulation phase, originally spatially situated
    at cell membrane midpoints.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCells
        Vector field caching the intracellular electric field.
    '''

    return VectorFieldCells(
        x=VectorCells(
            phase=phase, times_membranes_midpoint=phase.sim.efield_gj_x_time),
        y=VectorCells(
            phase=phase, times_membranes_midpoint=phase.sim.efield_gj_y_time),
    )


@type_check
def make_electric_extra(phase: SimPhase) -> VectorFieldCells:
    '''
    Vector field caching the electric field across all extracellular spaces for
    all time steps of the passed simulation phase, originally spatially situated
    at environmental grid space centres.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCells
        Vector field caching the extracellular electric field.

    Raises
    ----------
    BetseSimConfigException
        If this simulation has disabled extracellular spaces.
    '''

    # If extracellular spaces are disabled, raise an exception.
    phase.p.die_unless_ecm()

    # Create and return this field.
    return VectorFieldCells(
        x=VectorCells(
            phase=phase, times_grids_centre=phase.sim.efield_ecm_x_time),
        y=VectorCells(
            phase=phase, times_grids_centre=phase.sim.efield_ecm_y_time),
    )

# ....................{ MAKERS ~ microtubule               }....................
@type_check
def make_microtubule(phase: SimPhase) -> VectorFieldCells:
    '''
    Vector field caching all cellular microtubules for all time steps of the
    passed simulation phase, originally spatially situated at cell membrane
    midpoints.

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    Returns
    ----------
    VectorFieldCells
        Vector field caching all cellular microtubules.
    '''

    return VectorFieldCells(
        x=VectorCells(
            phase=phase, times_membranes_midpoint=phase.sim.mtubes_x_time),
        y=VectorCells(
            phase=phase, times_membranes_midpoint=phase.sim.mtubes_y_time),
    )
