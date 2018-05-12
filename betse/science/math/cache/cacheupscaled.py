#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level upscaled object subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from numpy import ndarray
from betse.science.export import expmath
from betse.science.math.cache.cacheabc import SimPhaseCacheABC
from betse.util.type.decorator.decmemo import property_cached

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

        return expmath.upscale_coordinates_tuple(
            self._phase.cells.xmin,
            self._phase.cells.xmax,
            self._phase.cells.ymin,
            self._phase.cells.ymax,
        )

    # ..................{ PROPERTIES ~ cells : vertices      }..................
    @property_cached
    def times_cells_vertices_coords(self) -> ndarray:
        '''
        Four-dimensional Numpy array of the upscaled coordinates of all cell
        membrane vertices for this cell cluster over all time steps under the
        assumption that these vertices are deformed over these time steps.

        See Also
        ----------
        :attrs:`betse.science.sim.Simulator.cell_verts_time`
            Further details.
        '''

        return expmath.upscale_coordinates(
            self._phase.sim.cell_verts_time)


    @property_cached
    def cells_vertices_coords(self) -> ndarray:
        '''
        Three-dimensional Numpy array of the upscaled coordinates of all cell
        membrane vertices for this cell cluster under the assumption that these
        vertices are *not* deformed at any time step and hence are fixed over
        all time steps.

        See Also
        ----------
        :attrs:`betse.science.cells.Cells.cell_verts`
            Further details.
        '''

        return expmath.upscale_coordinates(self._phase.cells.cell_verts)

    # ..................{ PROPERTIES ~ cells : centre        }..................
    @property_cached
    def cells_centre_x(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all cell
        centres for this cell cluster.
        '''

        return expmath.upscale_coordinates(
            self._phase.cells.cell_centres[:, 0])


    @property_cached
    def cells_centre_y(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled Y coordinates of all cell
        centres for this cell cluster.
        '''

        return expmath.upscale_coordinates(
            self._phase.cells.cell_centres[:, 1])

    # ..................{ PROPERTIES ~ grids : centre        }..................
    @property_cached
    def grids_centre_x(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all
        environmental grid space centres for this cell cluster.
        '''

        return expmath.upscale_coordinates(self._phase.cells.xypts[:, 0])


    @property_cached
    def grids_centre_y(self) -> ndarray:
        '''
        One-dimensional Numpy array of the upscaled X coordinates of all
        environmental grid space centres for this cell cluster.
        '''

        return expmath.upscale_coordinates(self._phase.cells.xypts[:, 1])
