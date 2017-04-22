#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level vector subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.simulate.cache.cacheabc import SimPhaseCacheABC
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheVectorCells(SimPhaseCacheABC):
    '''
    Simulation phase-specific vector subcache, persisting all previously
    constructed vectors for a single simulation phase.

    Attributes
    ----------
    layer : SimPhaseCacheLayerCellsVector
        Subcache of all vector-based layers constructed for this phase.
    '''

    # ..................{ INITIALIZORS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Avoid circular import dependencies.
        from betse.science.simulate.cache.cachelyrvec import (
            SimPhaseCacheLayerCellsVector)

        # Classify all subcaches imported above.
        self.layer = SimPhaseCacheLayerCellsVector(self._phase)

    # ..................{ PROPERTIES ~ currents              }..................
    @property_cached
    def voltage_membrane(self) -> VectorCellsCache:
        '''
        Vector cache of all upscaled **transmembrane voltages** (i.e., voltages
        across all gap junctions connecting intracellular membranes) over all
        sampled time steps of the current simulation phase, originally spatially
        situated at cell membrane midpoints.

        For readability of units in exported visuals (e.g., plots), voltages are
        cached upscaled from volts (V) to millivolts (mV).
        '''

        return VectorCellsCache(
            phase=self._phase,
            times_membranes_midpoint=expmath.upscale_cell_data(
                self._phase.sim.vm_time))

    # ..................{ PROPERTIES ~ deform                }..................
    #FIXME: Raise an exception unless deformations are enabled. To do so sanely,
    #we'll want to define a new @phase_property_cached decorator accepting an
    #optional "requirements" parameter, much like the existing @piperunner
    #decorator. (For now, simply ignore this for simplicity.)
    @property_cached
    def deform_total_magnitudes(self) -> VectorCellsCache:
        '''
        Vector cache of the upscaled magnitudes of all **total cellular
        deformations** (i.e., summations of all cellular deformations due to
        galvanotropic and osmotic pressure body forces) over all sampled time
        steps of the current simulation phase, originally spatially situated at
        cell membrane midpoints.
        '''

        # Two-dimensional Numpy arrays of the X and Y components of all total
        # cellular deformations over all time steps, mapped from cell centres to
        # cell membrane midpoints.
        times_membranes_x = (
            self._phase.sim.dx_cell_time[:, self.cells.mem_to_cells])
        times_membranes_y = (
            self._phase.sim.dy_cell_time[:, self.cells.mem_to_cells])

        # One-dimensional Numpy array of the magnitudes of all total
        # cellular deformations over all time steps, mapped from cell centres to
        # cell membrane midpoints.
        times_membranes_magnitudes = expmath.upscale_cells_coordinates(
            times_membranes_x * self._phase.cells.membrane_normal_unit_x +
            times_membranes_y * self._phase.cells.membrane_normal_unit_y)

        # Create, return, and cache this vector.
        return VectorCellsCache(
            phase=self._phase,
            times_membranes_midpoint=times_membranes_magnitudes)
