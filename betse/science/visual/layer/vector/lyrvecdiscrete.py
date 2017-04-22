#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading each individual cell in the cell cluster as a
discrete (i.e., individual, disseparate, discontiguous) surface.
'''

# ....................{ IMPORTS                            }....................
from abc import abstractproperty
from betse.science.visual.layer.vector.lyrvecabc import (
    LayerCellsVectorColorfulABC)
from betse.util.type.types import type_check, IterableTypes, SequenceOrNoneTypes
from numpy import ndarray

# ....................{ SUPERCLASSES                       }....................
class LayerCellsVectorDiscreteMembranesABC(LayerCellsVectorColorfulABC):
    '''
    Abstract base class of all layer subclasses plotting a single vector
    spatially situated at cell membrane vertices (e.g., transmembrane voltages)
    as a discontiguous Gouraud-shaded surface depicted by a polygonal mesh onto
    the cell cluster for one on more time steps.

    Such layers are somewhat more computationally expensive in both space and
    time than the average layer. Gouraud-shading the surface of each cell
    requires a triangulation mesh for that cell, which this layer internally
    computes via the Delaunay hull of that cell's non-triangular vertices.

    Attributes
    ----------
    _cell_tri_meshes : list
        List of :class:`matplotlib.collections.TriMesh` instances, each
        encapsulating the triangulation mesh for the cell in this cell cluster
        whose 0-based index is the same structure as that of the
        :attr:`betse.science.cells.Cells.cell_verts` array.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Default all instance attributes.
        self._cell_tri_meshes = None

    # ..................{ SUBCLASS                           }..................
    @abstractproperty
    def cells_vertices_coords(self) -> ndarray:
        '''
        Three-dimensional Numpy array of the upscaled coordinates of all cell
        membrane vertices for this cell cluster at the current time step.

        See Also
        ----------
        :attrs:`betse.science.cells.Cells.cell_verts`
            Further details.
        '''

        pass

    # ..................{ SUPERCLASS                         }..................
    @property
    def color_data(self) -> SequenceOrNoneTypes:
        return self._vector.times_membranes_vertex


    @type_check
    def _layer_first_color_mappables(self) -> IterableTypes:

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = self._vector.times_membranes_vertex[
            self._visual.time_step]

        # List of triangulation meshes created by iteration below.
        self._cell_tri_meshes = []

        # For the index and two-dimensional array of upscaled vertex X and Y
        # coordinates for each cell in this cluster...
        for cell_index, cell_vertices_coords in enumerate(
            self.cells_vertices_coords):
            # X coordinates of all vertices defining this cell's polygon.
            cell_vertices_x = cell_vertices_coords[:, 0]

            # Y coordinates of all vertices defining this cell's polygon.
            cell_vertices_y = cell_vertices_coords[:, 1]

            # Color values of all vertices of this cell, referred to as "C" in
            # both the documentation and implementation of the tripcolor()
            # function. Why "C"? Because you will believe.
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

# ....................{ SUBCLASSES                         }....................
class LayerCellsVectorDiscreteMembranesFixed(
    LayerCellsVectorDiscreteMembranesABC):
    '''
    Layer subclass plotting a single vector spatially situated at cell membrane
    vertices (e.g., transmembrane voltages) as a discontiguous Gouraud-shaded
    surface depicted by a polygonal mesh onto the cell cluster for one on more
    time steps under the assumption that these vertices are fixed over these
    time steps.

    This assumption dramatically improves the computational efficiency of this
    layer, permitting the initial set of cell triangulations computed for the
    first time step to be reused for all subsequent time steps. Naturally, this
    assumption breaks down for simulations enabling deformations.
    '''

    # ..................{ SUPERCLASS                         }..................
    @property
    def cells_vertices_coords(self) -> ndarray:
        return self._phase.cache.upscaled.cells_vertices_coords


    # For efficiency, this method simply reshades the triangulated mesh for each
    # cell previously computed by _layer_first_color_mappables().
    def _layer_next(self) -> None:

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex_data = self._vector.times_membranes_vertex[
            self._visual.time_step]

        # For the index and triangulation mesh for each cell...
        for cell_index, cell_tri_mesh in enumerate(self._cell_tri_meshes):
            # Color values of all membrane vertices of this cell.
            cell_membranes_vertex_data = membranes_vertex_data[
                self._phase.cells.cell_to_mems[cell_index]]

            # Gouraud-shade this triangulation mesh with these color values.
            cell_tri_mesh.set_array(cell_membranes_vertex_data)


class LayerCellsVectorDiscreteMembranesDeformed(
    LayerCellsVectorDiscreteMembranesABC):
    '''
    Layer subclass plotting a single vector spatially situated at cell membrane
    vertices (e.g., transmembrane voltages) as a discontiguous Gouraud-shaded
    surface depicted by a polygonal mesh onto the cell cluster for one on more
    time steps under the assumption that these vertices are deformed over these
    time steps.

    This assumption dramatically reduces the computational efficiency of this
    layer, requiring the initial set of cell triangulations computed for the
    first time step to be recomputed at each subsequent time steps. Naturally,
    this assumption is only required for simulations enabling deformations.
    '''

    # ..................{ SUPERCLASS                         }..................
    @property
    def cells_vertices_coords(self) -> ndarray:
        return self._phase.cache.upscaled.times_cells_vertices_coords[
            self._visual.time_step]


    def _layer_next(self) -> None:

        #FIXME: If this fails... we have no idea. Perhaps pilfer the current
        #matplotlib repository for a working remove method?  Note that the class
        #in question is "matplotlib.collections.TriMesh". Google shows no hits,
        #sadly. We may indeed need to entirely clear and reconstruct the current
        #figure axes (!) from the ground up, which would be... intense.

        # For the triangulation mesh previously computed for each cell...
        for cell_tri_mesh in self._cell_tri_meshes:
            # Remove this mesh from the current figure axes.
            cell_tri_mesh.remove()

        # Recompute and return triangulation meshes for all cells.
        return self._layer_first_color_mappables()
