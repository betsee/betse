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
            times_membranes_midpoint=expmath.upscale_units_milli(
                self._phase.sim.vm_time))
