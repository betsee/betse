#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all matplotlib-based layer subclasses spatially
plotting vector fields onto the cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.field.fieldcls import VectorFieldCells
from betse.science.visual.layer.layerabc import (
    LayerCellsABC, LayerCellsColorfulABC)
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class LayerCellsFieldColorlessABC(LayerCellsABC):
    '''
    Abstract base class of all classes spatially plotting vector fields of
    arbitrary data onto the cell cluster (e.g., current density, electric field)
    such that neither the components nor magnitudes of these fields are mappable
    as colors onto the colorbars of parent plots and animations.

    Attributes
    ----------
    _field : VectorFieldCells
        Cache of various vector fields of the same underlying data situated
        along different coordinate systems for all time steps to be animated.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, field: VectorFieldCells) -> None:
        '''
        Initialize this layer.

        Parameters
        ----------
        field : VectorFieldCells
            Cache of various vector fields of the same underlying data situated
            along different coordinate systems for all time steps to be
            animated.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify all passed parameters.
        self._field = field


#FIXME: Rename to "LayerCellsFieldColorfulABC".
class LayerCellsFieldColoredABC(LayerCellsColorfulABC):
    '''
    Abstract base class of all classes spatially plotting vector fields of
    arbitrary data onto the cell cluster (e.g., current density, electric field)
    such that the components and/or magnitudes of these fields are mappable as
    colors onto the colorbars of parent plots and animations.

    Attributes
    ----------
    _field : VectorFieldCells
        Cache of various vector fields of the same underlying data situated
        along different coordinate systems for all time steps to be animated.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, field: VectorFieldCells) -> None:
        '''
        Initialize this layer.

        Parameters
        ----------
        field : VectorFieldCells
            Cache of various vector fields of the same underlying data situated
            along different coordinate systems for all time steps to be
            animated.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify all passed parameters.
        self._field = field
