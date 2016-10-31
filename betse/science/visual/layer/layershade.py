#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsMappableABC
from betse.util.type.types import type_check, IterableTypes, SequenceTypes

# ....................{ CLASSES                            }....................
# class LayerCellsShadeContinuous(LayerCellsMappableABC):
#     '''
#     Layer plotting the entire cell cluster as a continuous Gouraud-shaded
#     surface represented as a polygonal mesh, interpolating the cell data for
#     each cell across the smooth spatial continuum inhabited by that cell.
#
#     Attributes
#     ----------
#     _cluster_tri_mesh : matplotlib.collections.TriMesh
#         Unstructured triangulation mesh grid interpolating the cell data for
#         each cell over the cell cluster as a contiuous whole.
#     '''
#
#     # ..................{ INITIALIZERS                       }..................
#     def __init__(self) -> None:
#         '''
#         Initialize this layer.
#         '''
#
#         # Initialize our superclass.
#         super().__init__()
#
#         # Default instance attributes.
#         self._cluster_tri_mesh = None
#
#     # ..................{ SUPERCLASS                         }..................
#     #FIXME: Implement us up.


class LayerCellsShadeDiscrete(LayerCellsMappableABC):
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
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._cell_tri_meshes = None

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:
        '''
        Layer the spatial distribution of a single cell-specific modelled
        variable (e.g., cell membrane voltage) for the first simulation time
        step onto the figure axes of the current plot or animation.

        Returns
        ----------
        IterableTypes
            Iterable of all mappables cached into the :attr:`_color_mappables`
            attribute by the :meth:`_layer_first` method.
        '''

        # Two-dimensional array of all cell data for the current time step.
        cells_vertex_data = self._get_cells_vertex_data()

        # Three-dimensional array of all upscaled cell vertex coordinates.
        # See "Cells.cell_verts" documentation for further details.
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

            # Average color values of all cell vertices, referred to as "C"
            # in both the documentation and implementation of the
            # tripcolor() function. Why "C"? Because you will believe.
            cell_vertex_data = cells_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shaded triangulation mesh for this cell, computed from
            # the Delaunay hull of the non-triangular vertices of this cell.
            cell_tri_mesh = self._visual.axes.tripcolor(
                # For equally obscure and uninteresting reasons, this
                # function requires these parameters be passed positionally.
                cell_vertices_x, cell_vertices_y, cell_vertex_data,

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
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the next time step and each cell of the current
        cluster onto the figure axes of the passed plot or animation as a
        discontiguous Gouraud-shaded surface represented as a polygonal mesh.

        For efficiency, this method simply reshades the triangulated mesh for
        each cell previously computed by the :meth:`_layer_first` method.
        '''

        # Two-dimensional array of all cell data for the current time step.
        cells_vertex_data = self._get_cells_vertex_data()

        # For the index and triangulation mesh for each cell...
        for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
            # Average color values of all cell vertices.
            cell_vertex_data = cells_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shade this triangulation mesh with these color values.
            cell_tri_mesh.set_array(cell_vertex_data)

    # ..................{ GETTERS                            }..................
    def _get_cells_vertex_data(self) -> SequenceTypes:
        '''
        Two-dimensional array of all cell data for the current time step,
        mapped from the midpoints onto the vertices of each cell membrane.
        '''

        return np.dot(
            self._visual.cell_data, self._visual.cells.matrixMap2Verts)
