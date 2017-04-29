#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
    # def prep(self, *args, **kwargs):
    #     super().prep(*args, **kwargs)
    #     import numpy as np
    #     _current_density_x_time_series = 100*np.asarray(self._phase.sim.I_tot_x_time)
    #     _current_density_y_time_series = 100*np.asarray(self._phase.sim.I_tot_y_time)
    #
    #     # Time series of all current density magnitudes (i.e., `Jmag_M`),
    #     # multiplying by 100 to obtain current density in units of uA/cm2.
    #     self._current_density_magnitude_time_series = np.sqrt(
    #         np.asarray(_current_density_x_time_series) ** 2 +
    #         np.asarray(_current_density_y_time_series) ** 2) + 1e-15
    #
    #     print('good[-1][0]: {}'.format(self._current_density_magnitude_time_series[-1][0]))
    #     print('badd[-1][0]: {}'.format(self._vector.times_grids_centre[-1][0]))


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

        # X and Y coordinates of the centers of all polygonol regions of the
        # Voronoi diagram defining this cell cluster.
        regions_centre_x = expmath.upscale_coordinates(
            self._phase.cells.voronoi_centres[:,0])
        regions_centre_y = expmath.upscale_coordinates(
            self._phase.cells.voronoi_centres[:,1])

        # One-dimensional array of all region-centred data for this time step.
        regions_centre_data = self._vector.times_regions_centre[
            self._visual.time_step]

        # Gouraud-shaded triangulation mesh for this cell cluster, computed from
        # the Delaunay hull of the non-triangular centers of these regions.
        self._cluster_tri_mesh = self._visual.axes.tripcolor(
            # Positional arguments. Thanks to internal flaws in the
            # matplotlib.tri.tripcolor() function parsing arguments passed
            # to the matplotlib.axes.tripcolor() method called here, the first
            # four arguments *MUST* be passed as positional arguments.
            regions_centre_x, regions_centre_y, regions_centre_data,

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

        # One-dimensional array of all region-centred data for this time step.
        regions_centre_data = self._vector.times_regions_centre[
            self._visual.time_step]

        # Gouraud-shade this triangulation mesh with these color values.
        self._cluster_tri_mesh.set_array(regions_centre_data)
