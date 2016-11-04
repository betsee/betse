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
from betse.util.type.types import type_check, IterableTypes  #, SequenceTypes
from numpy import ndarray

# ....................{ CLASSES                            }....................
#FIXME: Fix us up, please. This layer is effectively broken at the moment,
#plotting a spatially symmetric distribution even where the underlying data is
#asymmetric (in which case one would hope for some sort of distinct gradient).
#This layer is frequently leveraged elsewhere and hence fairly critical.

class LayerCellsShadeContinuous(LayerCellsMappableABC):
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
    def __init__(self) -> None:

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._cluster_tri_mesh = None

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # Upscaled X coordinates of the centers of all polygonol regions of the
        # Voronoi diagram defining this cell cluster.
        regions_x = visuals.upscale_cell_coordinates(
            self._visual.cells.voronoi_centres[:,0]),

        # Upscaled Y coordinates of the centers of all polygonol regions.
        regions_y = visuals.upscale_cell_coordinates(
            self._visual.cells.voronoi_centres[:,1]),

        # Average color values of all cell vertices, referred to as "C"
        # in both the documentation and implementation of the
        # tripcolor() function. Why "C"? Because you will believe.
        regions_data = self._get_regions_data()

        # Gouraud-shaded triangulation mesh for this cell cluster, computed from
        # the Delaunay hull of the non-triangular centers of these regions.
        self._cluster_tri_mesh = self._visual.axes.tripcolor(
            # For equally obscure and uninteresting reasons, this
            # function requires these parameters be passed positionally.
            regions_x, regions_y, regions_data,

            # All remaining parameters may be passed by keyword.
            shading='gouraud',
            cmap=self._visual.colormap,

            #FIXME: Minimum and maximum colors should absolutely be passed, but
            #presumably are *NOT* quite right at the moment. Shouldn't these
            #values be the minimum and maximum of the "regions_data" array?
            # vmin=self._visual.color_min,
            # vmax=self._visual.color_max,
        )

        # Map this triangulation mesh onto the figure colorbar.
        return self._cluster_tri_mesh


    def _layer_next(self) -> None:

        self._cluster_tri_mesh.set_array(self._get_regions_data())

    # ..................{ GETTERS                            }..................
    #FIXME: Generalize this into a new
    #AnimCellsMembranesData.regions_centre_data() property; then,
    #replace all calls to this method above by this property.

    def _get_regions_data(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all cell data for the current time step,
        mapped from the midpoints of all cell membranes onto the centres of all
        Voronoi regions.
        '''

        # Two-dimensional array of all cell data for the current time step.
        cells_centre_data = self._visual.cells_data

        #FIXME: Document us up.
        regions_centre_data = np.zeros(len(self._visual.cells.voronoi_centres))
        regions_centre_data[self._visual.cells.cell_to_grid] = cells_centre_data

        # Return this array.
        return regions_centre_data


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

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._cell_tri_meshes = None

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # Two-dimensional array of all cell data for the current time step.
        cells_membranes_vertex_data = self._get_cells_membranes_vertex_data()

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

            # Average color values of all cell vertices, referred to as "C"
            # in both the documentation and implementation of the
            # tripcolor() function. Why "C"? Because you will believe.
            cell_membranes_vertex_data = cells_membranes_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shaded triangulation mesh for this cell, computed from
            # the Delaunay hull of the non-triangular vertices of this cell.
            cell_tri_mesh = self._visual.axes.tripcolor(
                # For equally obscure and uninteresting reasons, this
                # function requires these parameters be passed positionally.
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

        # Two-dimensional array of all cell data for the current time step.
        cells_membranes_vertex_data = self._get_cells_membranes_vertex_data()

        # For the index and triangulation mesh for each cell...
        for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
            # Average color values of all cell vertices.
            cell_membranes_vertex_data = cells_membranes_vertex_data[
                self._visual.cells.cell_to_mems[cell_index]]

            # Gouraud-shade this triangulation mesh with these color values.
            cell_tri_mesh.set_array(cell_membranes_vertex_data)

    # ..................{ GETTERS                            }..................
    #FIXME: Generalize this into a new
    #AnimCellsMembranesData.cells_membranes_vertex_data() property; then,
    #replace all calls to this method above by this property.

    def _get_cells_membranes_vertex_data(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all cell data for the current time step,
        mapped from the midpoints onto the vertices of all cell membranes.
        '''

        return np.dot(
            self._visual.cells_membranes_data,

            #FIXME: Document this. Fairly intense, presumably.
            self._visual.cells.matrixMap2Verts)
