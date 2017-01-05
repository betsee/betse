#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
# import numpy as np
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsMappableArrayABC
from betse.util.type.types import type_check, IterableTypes, SequenceOrNoneTypes
# from numpy import ndarray

# ....................{ CLASSES                            }....................
#FIXME: Fix us up, please. This layer is effectively broken at the moment,
#plotting a spatially symmetric distribution even where the underlying data is
#asymmetric (in which case one would hope for some sort of distinct gradient).
#This layer is frequently leveraged elsewhere and hence fairly critical.

class LayerCellsShadeContinuous(LayerCellsMappableArrayABC):
    '''
    Layer plotting the entire cell cluster as a continuous Gouraud-shaded
    surface represented as a polygonal mesh, interpolating the cell data for
    each cell across the smooth spatial continuum inhabited by that cell.

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

        return self.times_regions_centre_data


    def _layer_first_color_mappables(self) -> IterableTypes:

        # Upscaled X coordinates of the centers of all polygonol regions of the
        # Voronoi diagram defining this cell cluster.
        regions_centre_x = visuals.upscale_cell_coordinates(
            self._visual.cells.voronoi_centres[:,0])

        # Upscaled Y coordinates of the centers of all polygonol regions.
        regions_centre_y = visuals.upscale_cell_coordinates(
            self._visual.cells.voronoi_centres[:,1])

        # One-dimensional array of all region centre data for this time step.
        regions_centre_data = (
            self.times_regions_centre_data[self._visual.time_step])

        # Gouraud-shaded triangulation mesh for this cell cluster, computed from
        # the Delaunay hull of the non-triangular centers of these regions.
        self._cluster_tri_mesh = self._visual.axes.tripcolor(
            # For equally obscure and uninteresting reasons, this
            # function requires these parameters to be passed positionally.
            regions_centre_x, regions_centre_y, regions_centre_data,

            # All remaining parameters may be passed by keyword.
            shading='gouraud',
            cmap=self._visual.colormap,
            vmin=self._visual.color_min,
            vmax=self._visual.color_max,
        )

        # Map this triangulation mesh onto the figure colorbar, returned as a
        # 1-tuple to comply with the superclass API.
        return (self._cluster_tri_mesh,)


    def _layer_next(self) -> None:

        # One-dimensional array of all region centre data for this time step.
        regions_centre_data = (
            self.times_regions_centre_data[self._visual.time_step])

        # Gouraud-shade this triangulation mesh with these color values.
        self._cluster_tri_mesh.set_array(regions_centre_data)


class LayerCellsShadeDiscrete(LayerCellsMappableArrayABC):
    '''
    Layer plotting each cell in the cell cluster as a discontiguous
    Gouraud-shaded surface represented as a polygonal mesh.

    This layer is somewhat more computationally expensive in both space and
    time than the average layer. Gouraud-shading the surface of each cell
    requires a triangulation mesh for that cell, which this layer internally
    computes via the Delaunay hull of that cell's non-triangular vertices.

    Attributes
    ----------
    _cell_tri_meshes : list
        List of :class:`matplotlib.collections.TriMesh` instances, each
        encapsulating the triangulation mesh for the cell in this cell cluster
        whose 0-based index is the same as that of the
        :attr:`betse.science.cells.Cells.cell_verts` array.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Default all instance attributes.
        self._cell_tri_meshes = None

    # ..................{ SUPERCLASS                         }..................
    @property
    def color_data(self) -> SequenceOrNoneTypes:

        return self.times_membranes_vertex_data


    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = (
            self.times_membranes_vertex_data[self._visual.time_step])

        # Three-dimensional array of all upscaled cell vertex coordinates. See
        # "Cells.cell_verts" documentation for further details.
        cells_vertices_coords = visuals.upscale_cell_coordinates(
            self._visual.cells.cell_verts)

        # List of triangulation meshes created by iteration below.
        self._cell_tri_meshes = []

        # For the index and two-dimensional array of vertex coordinates for
        # each cell in this cluster...
        for cell_index, cell_vertices_coords in enumerate(
            cells_vertices_coords):
            # X coordinates of all vertices defining this cell's polygon.
            cell_vertices_x = cell_vertices_coords[:, 0]

            # Y coordinates of all vertices defining this cell's polygon.
            cell_vertices_y = cell_vertices_coords[:, 1]

            # Average color values of all vertices of this cell, referred to as
            # "C" in both the documentation and implementation of the
            # tripcolor() function. Why "C"? Because you will believe.
            cell_membranes_vertex_data = membranes_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shaded triangulation mesh for this cell, computed from
            # the Delaunay hull of the non-triangular vertices of this cell.
            cell_tri_mesh = self._visual.axes.tripcolor(
                # For equally obscure and uninteresting reasons, this
                # function requires these parameters to be passed positionally.
                cell_vertices_x, cell_vertices_y, cell_membranes_vertex_data,

                # All remaining parameters may be passed by keyword.
                shading='gouraud',
                cmap=self._visual.colormap,
                vmin=self._visual.color_min,
                vmax=self._visual.color_max,
            )

            # Add this triangulation mesh to the cached set of such meshes.
            self._cell_tri_meshes.append(cell_tri_mesh)

        # Map these triangulation meshes onto the figure colorbar.
        return self._cell_tri_meshes


    def _layer_next(self) -> None:

        # For efficiency, this method simply reshades the triangulated mesh for
        # each cell previously computed by _layer_first_color_mappables().

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = (
            self.times_membranes_vertex_data[self._visual.time_step])

        # For the index and triangulation mesh for each cell...
        for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
            # Average color values of all vertices of this cell.
            cell_membranes_vertex_data = membranes_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shade this triangulation mesh with these color values.
            cell_tri_mesh.set_array(cell_membranes_vertex_data)

# ....................{ FACTORIES                          }....................
@type_check
def make(p: 'betse.science.parameters.Parameters', *args, **kwargs) -> (
    LayerCellsMappableArrayABC):
    '''
    Layer plotting the cell cluster as a Gouraud-shaded surface in either a
    contiguous or discontiguous manner according to the passed configuration.

    Parameters
    ----------
    p : Parameters
        Current simulation configuration.

    All remaining parameters are passed to either the
    :class:`LayerCellsShadeDiscrete` or :class:`LayerCellsShadeContinuous`
    constructor as is.
    '''

    # Type of layer to be created.
    layer_type = (
        LayerCellsShadeDiscrete if p.showCells else
        LayerCellsShadeContinuous)

    # Create and return an instance of this type, passed the passed parameters.
    return layer_type(*args, **kwargs)
