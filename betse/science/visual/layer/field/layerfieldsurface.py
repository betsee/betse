#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying arbitrary extracellular data as
interpolated surfaces (e.g., contour maps, heat maps) onto the cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.visual.layer.field.layerfieldabc import (
    LayerCellsFieldColoredABC)
from betse.util.type.types import type_check, IterableTypes

# ....................{ SUBCLASSES                         }....................
class LayerCellsFieldSurface(LayerCellsFieldColoredABC):
    '''
    Layer subclass plotting the magnitudes of a single vector field (e.g.,
    electric field) interpolated as a continuous surface onto the centres of all
    environmental grid spaces for one on more simulation time steps.

    Attributes
    ----------
    _surface_plot : matplotlib.image.AxesImage
        Surface plot of all magnitudes of this field previously plotted for the
        prior time step if any or `None` otherwise.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass.
        super().__init__(*args, **kwargs)

        # Default all remaining instance variables.
        self._surface_plot = None

    # ..................{ SUPERCLASS                         }..................
    def _layer_first_color_mappables(self) -> IterableTypes:
        '''
        Layer the interpolated magnitudes of this vector field onto the first
        simulation time step onto the figure axes of the current plot or
        animation, returning the mappable or mappables with which to map color
        data onto the colorbar.
        '''

        # Two-dimensional Numpy array of all field magnitudes for this time step
        # spatially situated at environmental grid space centres.
        field_magnitudes = self._field.times_grids_centre.magnitudes[
            self._visual.time_step]

        # 4-tuple of the upscaled minimum and maximum X and Y coordinates for
        # this environmental grid.
        cells_extent = expmath.upscale_cells_coordinates(
            self._visual.cells.xmin,
            self._visual.cells.xmax,
            self._visual.cells.ymin,
            self._visual.cells.ymax,
        )

        # Surface image plot of these magnitudes plotted for this time step. See
        # the matplotlib.axes.imshow() docstring for further details.
        self._surface_plot = self._visual.axes.imshow(
            # Two- or three-dimensional Numpy array to be interpolated.
            X=field_magnitudes,

            # Colormap converting input data values into output color values.
            cmap=self._visual.colormap,

            # X and Y coordinates for the boundaries of this environmental grid.
            extent=cells_extent,

            # Machine-readable string specifying the type of interpolation
            # to perform, established by the current configuration.
            interpolate=self._p.interp_type,

            # The [0, 0] index of this two- or three-dimensional Numpy array
            # resides at the lower-left corner of this figure's axes.
            origin='lower',
        )

        # Map this surface image plot onto the figure colorbar, returned as a
        # 1-tuple to comply with the superclass API.
        return (self._surface_plot,)


    def _layer_next(self) -> None:
        '''
        Layer the most significant X and Y components of this vector field for
        the next time step onto the figure axes of the current plot or
        animation.
        '''

        # Two-dimensional Numpy array of all field magnitudes for this time step
        # spatially situated at environmental grid space centres.
        field_magnitudes = self._field.times_grids_centre.magnitudes[
            self._visual.time_step]

        # Replace all normalized X and Y components previously plotted for the
        # prior time step by these components.
        self._surface_plot.set_data(field_magnitudes)
