#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-based plotter of shaded cell clusters.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.plot.plot import upscale_data
from betse.science.plot.plotter.plotterabc import PlotterCellsABC
from betse.util.type.types import type_check

# ....................{ BASE                               }....................
class PlotterCellsGouraudShaded(PlotterCellsABC):
    '''
    Plotter subclass spatially plotting each cell in this cell cluster as a
    discontiguous Gouraud-shaded surface represented as a polygonal mesh.

    This plotter is somewhat more computationally expensive in both space and
    time than the average plotter. Gouraud-shading the surface of each cell
    requires a triangulation mesh for that cell, which this plotter internally
    computes via the Delaunay hull of that cell's non-triangular vertices.

    Attributes
    ----------
    _cell_tri_meshes : set
        Set of `matplotlib.collections.TriMesh` instances, each encapsulating
        the triangulation mesh for an arbitrary cell in this cell cluster.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this plotter.
        '''

        # Initialize our superclass.
        super().__init__()

        # Sanitize all instance attributes.
        self._cell_tri_meshes = None

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def plot(self, plot: 'betse.science.plot.PlotCellsABC') -> None:
        '''
        Plot a single modelled variable (e.g., membrane voltage) for each cell
        of this cell cluster as a discontiguous Gouraud-shaded surface
        represented as a polygonal mesh onto the figure axes of the passed
        parent plot or animation for the current simulation time step.

        Parameters
        ----------
        plot : PlotCellsABC
            Parent plot or animation instance to plot onto.
        '''

        # Two-dimensional array of all cell data for the current time step,
        # mapped from the midpoints onto the vertices of each cell membrane.
        cells_vertex_data = np.dot(plot.cell_data, plot.cells.matrixMap2Verts)

        # If the cell cluster has yet to be triangulated, this *MUST* be the
        # first call to this method. In this case, create this triangulation.
        if self._cell_tri_meshes is None:
            # Set of triangulation meshes created by iteration below.
            self._cell_tri_meshes = set()

            # Three-dimensional array of all upscaled cell vertex coordinates.
            # See "Cells.cell_verts" documentation for further details.
            cells_vertices_coords = upscale_data(plot.cells.cell_verts)

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
                    plot.cells.cell_to_mems[cell_index]]

                # Gouraud-shaded triangulation mesh for this cell, computed from
                # the Delaunay hull of the non-triangular vertices of this cell.
                cell_tri_mesh = plot.axes.tripcolor(
                    # For equally obscure and uninteresting reasons, this
                    # function requires these parameters be passed positionally.
                    cell_vertices_x, cell_vertices_y, cell_vertex_data,

                    # All remaining parameters may be passed by keyword.
                    shading='gouraud',
                    cmap=plot.colormap,
                    vmin=plot.color_min,
                    vmax=plot.color_max,
                )

                # Add this triangulation mesh to the cached set of such meshes.
                self._cell_tri_meshes.add(cell_tri_mesh)
        # Else, the cell cluster has already been triangulated. In this case,
        # simply reshade the triangulated mesh for each cell.
        else:
            # For the index and triangulation mesh for each cell...
            for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
                # Average color values of all cell vertices.
                cell_vertex_data = cells_vertex_data[
                    plot.cells.cell_to_mems[cell_index]]

                # Gouraud-shade this triangulation mesh with these color values.
                cell_tri_mesh.set_array(cell_vertex_data)
