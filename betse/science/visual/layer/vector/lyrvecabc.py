#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all matplotlib-based layer subclasses spatially
plotting vectors onto the cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.visual.layer.lyrabc import LayerCellsColorfulABC
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class LayerCellsVectorColorfulABC(LayerCellsColorfulABC):
    '''
    Abstract base class of all classes spatially plotting vectors of arbitrary
    data onto the cell cluster (e.g., cell voltage) such that the elements of
    these vectors are mappable as colors onto the colorbars of parent plots and
    animations.

    Attributes
    ----------
    _vector : VectorCellsCache
        Cache of various vectors of the same underlying data situated along
        different coordinate systems for all time steps to be animated.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, vector: VectorCellsCache) -> None:
        '''
        Initialize this layer.

        Parameters
        ----------
        vector : VectorCellsCache
            Cache of various vectors of the same underlying data situated along
            different coordinate systems for all time steps to be animated.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify all passed parameters.
        self._vector = vector
