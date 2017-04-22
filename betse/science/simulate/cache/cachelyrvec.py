#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level vector subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.simulate.cache.cacheabc import SimPhaseCacheABC
from betse.science.visual.layer.lyrabc import LayerCellsABC
from betse.science.visual.layer.vector.lyrvecdiscrete import (
    LayerCellsVectorDiscreteMembranesFixed)
from betse.science.visual.layer.vector.lyrvecsmooth import (
    LayerCellsVectorSmoothRegions)
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheLayerCellsVector(SimPhaseCacheABC):
    '''
    Simulation phase-specific **vector layer** (i.e., vector-based layer
    layering a one-dimensional Numpy array onto exported visuals) subcache,
    persisting all previously constructed vector-based layers for a single
    simulation phase.
    '''

    # ..................{ PROPERTIES ~ currents              }..................
    @property_cached
    def voltage_membrane(self) -> LayerCellsABC:
        '''
        Vector layer layering all transmembrane voltages (Vmem) for the cell
        cluster over all time steps onto caller-provided visuals.
        '''

        # Type of layer to be created, plotting the cell cluster as a
        # Gouraud-shaded surface in either a contiguous or discontiguous manner
        # according to the phase configuration.
        if self._phase.p.showCells:
            layer_type = LayerCellsVectorDiscreteMembranesFixed
        else:
            layer_type = LayerCellsVectorSmoothRegions

        # Return the desired layer.
        return layer_type(vector=self._phase.cache.vector.voltage_membrane)
