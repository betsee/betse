#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading each individual cell in the cell cluster as a
discontiguous surface.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.visual.layer.vector.lyrvecabc import (
    LayerCellsVectorColorfulABC)
from betse.util.type.types import type_check, IterableTypes, SequenceOrNoneTypes

# ....................{ SUBCLASSES                         }....................
class LayerCellsVectorAbruptMembranes(LayerCellsVectorColorfulABC):
    '''
    Layer subclass plotting a single vector spatially situated at cell membrane
    vertices (e.g., transmembrane voltages) as a discontiguous Gouraud-shaded
    surface depicted by a polygonal mesh onto the cell cluster for one on more
    simulation steps.

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

        return self._vector.times_membranes_vertex


    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = self._vector.times_membranes_vertex[
            self._visual.time_step]

        # Three-dimensional array of all upscaled cell vertex coordinates. See
        # "Cells.cell_verts" documentation for further details.
        cells_vertices_coords = expmath.upscale_cell_coordinates(
            self._phase.cells.cell_verts)

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
                self._phase.cells.cell_to_mems[cell_index]]

            # Gouraud-shaded triangulation mesh for this cell, computed from
            # the Delaunay hull of the non-triangular vertices of this cell.
            cell_tri_mesh = self._visual.axes.tripcolor(
                # Positional arguments. Thanks to internal flaws in the
                # matplotlib.tri.tripcolor() function parsing arguments passed
                # to the matplotlib.axes.tripcolor() method called here, the
                # first four arguments *MUST* be passed as positional arguments.
                cell_vertices_x, cell_vertices_y, cell_membranes_vertex_data,

                # Keyword arguments. All remaining arguments *MUST* be passed as
                # keyword arguments.
                shading='gouraud',
                vmin=self._visual.color_min,
                vmax=self._visual.color_max,

                # Colormap converting input values into output color values.
                cmap=self._visual.colormap,
            )

            # Add this triangulation mesh to the cached set of such meshes.
            self._cell_tri_meshes.append(cell_tri_mesh)

        # Map these triangulation meshes onto the figure colorbar.
        return self._cell_tri_meshes


    # For efficiency, this method simply reshades the triangulated mesh for each
    # cell previously computed by _layer_first_color_mappables().
    def _layer_next(self) -> None:

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = self._vector.times_membranes_vertex[
            self._visual.time_step]

        # For the index and triangulation mesh for each cell...
        for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
            # Average color values of all vertices of this cell.
            cell_membranes_vertex_data = membranes_vertex_data[
                self._phase.cells.cell_to_mems[cell_index]]

            # Gouraud-shade this triangulation mesh with these color values.
            cell_tri_mesh.set_array(cell_membranes_vertex_data)


