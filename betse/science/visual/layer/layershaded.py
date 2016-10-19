#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading the current cell cluster.
'''

#FIXME: Rename this submodule to "layershade".

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.types import type_check, IterableTypes

# ....................{ BASE                               }....................
#FIXME: Rename to "LayerCellsShadeDiscrete", in anticipation of a
#new "LayerCellsShadeContinuum" subclass.
class LayerCellsGouraudShaded(LayerCellsABC):
    '''
    Layer subclass spatially plotting each cell in the current cell cluster as a
    discontiguous Gouraud-shaded surface represented as a polygonal mesh.

    This layer is somewhat more computationally expensive in both space and
    time than the average layer. Gouraud-shading the surface of each cell
    requires a triangulation mesh for that cell, which this layer internally
    computes via the Delaunay hull of that cell's non-triangular vertices.

    Attributes
    ----------
    _cell_tri_meshes : list
        List of `matplotlib.collections.TriMesh` instances, each encapsulating
        the triangulation mesh for the cell in this cell cluster whose 0-based
        index is the same as that of the
        :attr:`betse.science.cells.Cells.cell_verts` array.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Sanitize all instance attributes.
        self._cell_tri_meshes = None

    # ..................{ PROPERTIES                         }..................
    #FIXME: This property is somewhat but *NOT* perfectly general. It doesn't
    #particularly apply to the perfecetly general-purpose "LayerCellsABC"
    #class (e.g., because of text-only plotters, for example). Perhaps a new
    #"PlotterCellsMappableABC" class should be constructed, in which case this
    #property should be shifted there. The class hiercharcy in that case would
    #resemble:
    #
    #    LayerCellsABC
    #           ^
    #           |
    #    PlotterCellsMappableABC
    #           ^
    #           |
    #    LayerCellsGouraudShaded

    @property
    def color_mappables(self) -> IterableTypes:
        '''
        Iterable of all :class:`matplotlib.cm.ScalarMappable` instances plotted
        by this layer.

        Note that, due to Matplotlib constraints, only the first mappable in
        this iterable defines the range of colors and hence colorbar of the
        parent plot or animation.
        '''

        return self._cell_tri_meshes

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def layer(
        self, visual: 'betse.science.visual.visualabc.VisualCellsABC') -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the current time step and each cell of the current
        cluster onto the figure axes of the passed plot or animation as a
        discontiguous Gouraud-shaded surface represented as a polygonal mesh.

        Parameters
        ----------
        visual : VisualCellsABC
            Plot or animation to layer onto.
        '''

        #FIXME: Do we compute this elsewhere in the codebase? If so, globally
        #cache this somewhere for reuse.

        # Two-dimensional array of all cell data for the current time step,
        # mapped from the midpoints onto the vertices of each cell membrane.
        cells_vertex_data = np.dot(
            visual.cell_data, visual.cells.matrixMap2Verts)

        # If the cell cluster has yet to be triangulated, this *MUST* be the
        # first call to this method. In this case, create this triangulation.
        if self._cell_tri_meshes is None:
            # List of triangulation meshes created by iteration below.
            self._cell_tri_meshes = []

            #FIXME: Do we compute this elsewhere in the codebase? If so,
            #globally cache this somewhere for reuse.

            # Three-dimensional array of all upscaled cell vertex coordinates.
            # See "Cells.cell_verts" documentation for further details.
            cells_vertices_coords = visuals.upscale_cell_coordinates(
                visual.cells.cell_verts)

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
                    visual.cells.cell_to_mems[cell_index]]

                # Gouraud-shaded triangulation mesh for this cell, computed from
                # the Delaunay hull of the non-triangular vertices of this cell.
                cell_tri_mesh = visual.axes.tripcolor(
                    # For equally obscure and uninteresting reasons, this
                    # function requires these parameters be passed positionally.
                    cell_vertices_x, cell_vertices_y, cell_vertex_data,

                    # All remaining parameters may be passed by keyword.
                    shading='gouraud',
                    cmap=visual.colormap,
                    vmin=visual.color_min,
                    vmax=visual.color_max,
                )
                # cell_tri_mesh.set_clim(visual.color_min, visual.color_max)

                # Add this triangulation mesh to the cached set of such meshes.
                self._cell_tri_meshes.append(cell_tri_mesh)
        # Else, the cell cluster has already been triangulated. In this case,
        # simply reshade the triangulated mesh for each cell.
        else:
            # For the index and triangulation mesh for each cell...
            for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
                # Average color values of all cell vertices.
                cell_vertex_data = cells_vertex_data[
                    visual.cells.cell_to_mems[cell_index]]

                # Gouraud-shade this triangulation mesh with these color values.
                cell_tri_mesh.set_array(cell_vertex_data)
