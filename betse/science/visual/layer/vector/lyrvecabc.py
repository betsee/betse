#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all matplotlib-based layer subclasses spatially
plotting vectors onto the cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.layer.lyrabc import LayerCellsColorfulABC
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class LayerCellsVectorColorfulABC(LayerCellsColorfulABC):
    '''
    Abstract base class of all classes spatially plotting vectors of arbitrary
    one-dimensional data (e.g., transmembrane voltage) onto the cell cluster
    such that the elements of these vectors are mappable as colors onto the
    colorbars of parent visuals.

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

# ....................{ MAKERS                             }....................
#FIXME: Rename to layer_vector().
@type_check
def make_layer(phase: SimPhase, *args, **kwargs) -> LayerCellsVectorColorfulABC:
    '''
    Vector Layer suitable for spatially plotting the passed vector of arbitrary
    one-dimensional data (e.g., transmembrane voltage) for the passed simulation
    phase onto the cell cluster such that the elements of this vector are
    mappable as colors onto the colorbars of parent visuals.

    The type of layer instantiated and returned by this factory function
    conditionally depends on the configuration for this phase. (See this
    function's implementation for details.)

    Parameters
    ----------
    phase : SimPhase
        Current simulation phase.

    All remaining arguments are passed as is to the ``__init__`` method of this
    layer's class.

    Returns
    ----------
    LayerCellsVectorColorfulABC
        Layer suitable for plotting the passed vector.
    '''

    # Avoid circular import dependencies.
    from betse.science.visual.layer.vector.lyrvecdiscrete import (
        LayerCellsVectorDiscreteMembranesFixed,
        LayerCellsVectorDiscreteMembranesDeformed,
    )
    from betse.science.visual.layer.vector.lyrvecsmooth import (
        LayerCellsVectorSmoothRegions)

    # Type of layer to be created.
    layer_type = None

    # If this simulation configuration requests that each cell in the cell
    # cluster be plotted...
    if phase.p.showCells:
        # If this phase enables deformations, dynamically plot these cells
        # in a manner redrawing cell boundaries each time step.
        if phase.p.deformation:
            layer_type = LayerCellsVectorDiscreteMembranesDeformed
        # Else, this phase does *NOT* enable deformations. In this case,
        # statically plot these cells in a manner reusing the cell
        # boundaries plotted for the first time step on subsequent time
        # steps.
        else:
            layer_type = LayerCellsVectorDiscreteMembranesFixed
    # Else, this simulation configuration requests that the cell cluster be
    # plotted as a smooth continuum. In this case, do so.
    else:
        layer_type = LayerCellsVectorSmoothRegions

    # Create and return this layer.
    return layer_type(*args, **kwargs)
