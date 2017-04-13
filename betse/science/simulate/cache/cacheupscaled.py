#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level upscaled object subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.simulate.cache.cacheabc import SimPhaseCacheABC
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check
from numpy import ndarray

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheUpscaled(SimPhaseCacheABC):
    '''
    Simulation phase-specific upscaled object subcache, persisting all
    previously constructed objects whose X and Y coordinates have been upscaled
    for use with :mod:`matplotlib` for a single simulation phase.

    Design
    ----------
    All objects persisted by this subcache are required only sporadically by
    one or more isolated features (e.g., plots, animations). Hence, these
    objects are *not* suitable for ownership by the general-purpose
    :class:`Cells` class.
    '''

    # ..................{ PROPERTIES                         }..................
    @property_cached
    def extent(self) -> tuple:
        '''
        4-tuple of the upscaled minimum and maximum X and Y coordinates of the
        environmental grid for this cell cluster.

        Specifically, this tuple consists of (in order):

        #. This grid's upscaled minimum X coordinate.
        #. This grid's upscaled maximum X coordinate.
        #. This grid's upscaled minimum Y coordinate.
        #. This grid's upscaled maximum Y coordinate.
        '''

        return expmath.upscale_cells_coordinates(
            self._phase.cells.xmin,
            self._phase.cells.xmax,
            self._phase.cells.ymin,
            self._phase.cells.ymax,
        )

    # ..................{ PROPERTIES ~ cells centre          }..................
    @property_cached
    def cells_centre_x(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all cell
        centres for this cell cluster.
        '''

        return expmath.upscale_cells_coordinates(
            self._phase.cells.cell_centres[:,0])


    @property_cached
    def cells_centre_y(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled Y coordinates of all cell
        centres for this cell cluster.
        '''

        return expmath.upscale_cells_coordinates(
            self._phase.cells.cell_centres[:,1])

    # ..................{ PROPERTIES ~ grids centre          }..................
    @property_cached
    def grids_centre_x(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all
        environmental grid space centres for this cell cluster.
        '''

        return expmath.upscale_cells_coordinates(self._phase.cells.xypts[:,0])


    @property_cached
    def grids_centre_y(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all
        environmental grid space centres for this cell cluster.
        '''

        return expmath.upscale_cells_coordinates(self._phase.cells.xypts[:,1])
