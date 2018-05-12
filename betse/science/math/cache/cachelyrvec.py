#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level vector subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.cache.cacheabc import SimPhaseCacheABC
from betse.science.visual.layer.vector import lyrvecabc
from betse.science.visual.layer.vector.lyrvecabc import (
    LayerCellsVectorColorfulABC)
from betse.util.type.decorator.decmemo import property_cached

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
    def voltage_membrane(self) -> LayerCellsVectorColorfulABC:
        '''
        Vector layer layering all transmembrane voltages (Vmem) for the cell
        cluster over all time steps onto caller-provided visuals.
        '''

        return lyrvecabc.make_layer(
            phase=self._phase, vector=self._phase.cache.vector.voltage_membrane)
