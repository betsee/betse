#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation phase cache** (e.g., container persisting previously
constructed large-scale objects for a simulation phase) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.vectorcls import VectorCells
from betse.science.math.vector.field.fieldcls import VectorFieldCells
from betse.science.simulate.cache.cacheabc import SimPhaseSubcacheABC
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
#FIXME: Refactor all "betse.science.math.vector.vectormake" and
#""betse.science.math.vector.field.fieldmake" functions into cached properties
#of this class.
class SimPhaseSubcacheVector(SimPhaseSubcacheABC):
    '''
    Lower-level simulation phase vector subcache, persisting all previously
    constructed vectors and vector fields for a single simulation phase.
    '''

    # ..................{ PROPERTIES                         }..................
    @property_cached
    def voltage_membrane_polarity(self) -> VectorFieldCells:
        '''
        Vector field cache of all transmembrane voltage (Vmem) polarities over
        all time steps of the current simulation phase, originally spatially
        situated at cell membrane midpoints.
        '''

        #FIXME: Document us up.
        #FIXME: Is "sim.vm_ave_time" a list? If so, we'll need to Numpify here.
        polm = (
            self._phase.sim.vm -
            self._phase.sim.vm_ave_time[:,self._phase.cells.mem_to_cells])

        #FIXME: What does multiplication by "polm" do here?
        # X and Y coordinates of all normal unit vectors orthogonal to all
        # tangent unit vectors of all cell membranes.
        polx = polm*self._phase.cells.mem_vects_flat[:,2]
        poly = polm*self._phase.cells.mem_vects_flat[:,3]

        #FIXME: Document us up.
        pcx = self._phase.map_membranes_midpoint_to_cells_centre(
            polx*self._phase.cells.mem_sa) / self._phase.cells.cell_sa
        pcy = self._phase.map_membranes_midpoint_to_cells_centre(
            poly*self._phase.cells.mem_sa) / self._phase.cells.cell_sa

        # Create, return, and cache this vector field.
        return VectorFieldCells(
            x=VectorCells(phase=self._phase, times_cells_centre=pcx),
            y=VectorCells(phase=self._phase, times_cells_centre=pcy))
