#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading the cell cluster as a continuous surface.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.visual.layer.vector.lyrvecabc import (
    LayerCellsVectorColorfulABC)
from betse.util.type.types import type_check, IterableTypes, SequenceOrNoneTypes

# ....................{ SUBCLASSES                         }....................
class LayerCellsVectorSmoothGrids(LayerCellsVectorColorfulABC):
    '''
    Layer subclass plotting a single vector spatially situated at environmental
    grid space centres (e.g., extracellular electric field magnitudes) as a
    continuous surface depicted by a polygonal mesh onto the cell cluster for
    one on more simulation steps.

    Attributes
    ----------
    _surface_image : matplotlib.image.AxesImage
        Surface image plot of this vector for the current time step.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Default all instance variables.
        self._surface_image = None

    # ..................{ SUPERCLASS                         }..................
    @property
    def color_data(self) -> SequenceOrNoneTypes:
        return self._vector.times_grids_centre


    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # Surface image plot of these magnitudes plotted for this time step. See
        # the matplotlib.axes.imshow() docstring for further details.
        self._surface_image = self._visual.axes.imshow(
            # Two-dimensional array of all grid data for this time step,
            # spatially situated at environmental grid space centres.
            X=self._vector.times_grids_centre[self._visual.time_step],
            # self._current_density_magnitude_time_series[self._visual.time_step],

            # Colormap converting input data values into output color values.
            cmap=self._visual.colormap,

            # X and Y coordinates of the boundaries of this environmental grid.
            extent=self._phase.cache.upscaled.extent,

            # The [0, 0] index of this two- or three-dimensional Numpy array
            # resides at the lower-left corner of this figure's axes.
            origin='lower',

            # Z-order of this plot with respect to other artists.
            zorder=self._zorder,
        )
        # self._surface_image.set_clim(0, 152)

        # Map this surface image plot onto the figure colorbar, returned as a
        # 1-tuple to comply with the superclass API.
        return (self._surface_image,)


    def _layer_next(self) -> None:

        # print('Here!')
        # self._visual._rescale_color_mappables()

        # Replace all obsoleted grid data plotted for the prior time step by the
        # updated grid data for this time step.
        self._surface_image.set_data(
            # self._current_density_magnitude_time_series[-1])
            # self._current_density_magnitude_time_series[self._visual.time_step])
            self._vector.times_grids_centre[self._visual.time_step])


class LayerCellsVectorSmoothRegions(LayerCellsVectorColorfulABC):
    '''
    Layer subclass plotting a single vector spatially situated at Voronoi region
    centres (e.g., transmembrane voltage averages) as a continuous
    Gouraud-shaded surface depicted by a polygonal mesh onto the cell cluster
    for one on more simulation steps, interpolating the cell data for each cell
    across the smooth spatial continuum inhabited by that cell.

    Attributes
    ----------
    _cluster_tri_mesh : matplotlib.collections.TriMesh
        Unstructured triangulation mesh grid interpolating the cell data for
        each cell over the cell cluster as a contiuous whole.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Default all instance attributes.
        self._cluster_tri_mesh = None

    # ..................{ SUPERCLASS                         }..................
    @property
    def color_data(self) -> SequenceOrNoneTypes:

        return self._vector.times_regions_centre


    def _layer_first_color_mappables(self) -> IterableTypes:

        # X and Y coordinates of all cell membrane midpoints.
        membranes_midpoint_x = expmath.upscale_coordinates(
            self._phase.cells.mem_mids_flat[:,0])
        membranes_midpoint_y = expmath.upscale_coordinates(
            self._phase.cells.mem_mids_flat[:,1])

        # Membrane midpoint-centred data for this time step.
        membranes_midpoint_data = self._vector.times_membranes_midpoint[
            self._visual.time_step]

        # Gouraud-shaded triangulation mesh for this cell cluster, computed from
        # the Delaunay hull of the non-triangular centers of these regions.
        self._cluster_tri_mesh = self._visual.axes.tripcolor(
            # Positional arguments. Thanks to internal flaws in the
            # matplotlib.tri.tripcolor() function parsing arguments passed
            # to the matplotlib.axes.tripcolor() method called here, the first
            # four arguments *MUST* be passed as positional arguments.
            membranes_midpoint_x, membranes_midpoint_y, membranes_midpoint_data,

            # Keyword arguments. All remaining arguments *MUST* be passed as
            # keyword arguments.
            shading='gouraud',
            vmin=self._visual.color_min,
            vmax=self._visual.color_max,

            # Colormap converting input data values into output color values.
            cmap=self._visual.colormap,

            # Z-order of this mesh with respect to other artists.
            zorder=self._zorder,
        )

        # Map this triangulation mesh onto the figure colorbar, returned as a
        # 1-tuple to comply with the superclass API.
        return (self._cluster_tri_mesh,)


    def _layer_next(self) -> None:

        # Membrane midpoint-centred data for this time step.
        membranes_midpoint_data = self._vector.times_membranes_midpoint[
            self._visual.time_step]

        # Gouraud-shade this triangulation mesh with these color values.
        self._cluster_tri_mesh.set_array(membranes_midpoint_data)