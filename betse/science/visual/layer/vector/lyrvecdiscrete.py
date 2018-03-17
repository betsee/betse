#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially shading each individual cell in the cell cluster as a
discrete (i.e., individual, disseparate, discontiguous) surface.
'''

# ....................{ IMPORTS                            }....................
from betse.science.visual.layer.vector.lyrvecabc import (
    LayerCellsVectorColorfulABC)
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.types import type_check, IterableTypes, SequenceOrNoneTypes
from matplotlib.collections import TriMesh
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
        # self._visual.axes.set_axis_bgcolor('black')

        # One-dimensional array of all membrane vertex data for this time step.
        membranes_vertex = self._vector.times_membranes_vertex[
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

            # Indices of all membranes for this cell.
            cell_membranes_index = self._phase.cells.cell_to_mems[cell_index]

            # Color values of all vertices of this cell, referred to as "C" in
            # both the documentation and implementation of the tripcolor()
            # function. Why "C"? Because you will believe.
            cell_membranes_vertex = membranes_vertex[cell_membranes_index]

            # Gouraud-shaded triangulation mesh for this cell, computed from
            # the Delaunay hull of the non-triangular vertices of this cell.
            cell_tri_mesh = self._visual.axes.tripcolor(
                # Positional arguments. Thanks to internal flaws in the
                # matplotlib.tri.tripcolor() function parsing arguments passed
                # to the matplotlib.axes.tripcolor() method called here, the
                # first four arguments *MUST* be passed as positional arguments.
                cell_vertices_x, cell_vertices_y, cell_membranes_vertex,

                # Keyword arguments. All remaining arguments *MUST* be passed as
                # keyword arguments.
                shading='gouraud',
                vmin=self._visual.color_min,
                vmax=self._visual.color_max,

                # Colormap converting input values into output color values.
                cmap=self._visual.colormap,

                # Z-order of this mesh with respect to other artists.
                zorder=self._zorder,
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


    # Of necessity, this method recomputes the triangulated mesh for each
    # cell by recalling _layer_first_color_mappables() each time step. To ensure
    # that these meshes are properly associated with the existing colorbar, the
    # _layer_next_color_mappables() rather than _layer_next() method is defined.
    def _layer_next_color_mappables(self) -> None:

        #FIXME: Submit a matplotlib issue concerning this, please.

        # Manually remove all triangulation meshes plotted for the prior time
        # step by iterating over all matplotlib collection objects and
        # preserving all such objects that are *NOT* such meshes. Doing so also
        # removes all triangulation meshes of other visuals already plotted for
        # this time step and is hence non-ideal. But no alternatives exist.
        #
        # Actually, that's not quite accurate. The following alternative exists,
        # but silently reduces to a noop for unknown reasons:
        #
        #     # For each triangulation mesh plotted for the prior time step...
        #     for cell_tri_mesh in self._cell_tri_meshes:
        #         # Remove this mesh from the current figure axes. Sadly, doing
        #         # so via this method call appears to silently reduce to a noop
        #         # rather than raising the expected "NotImplementedError"
        #         # exception. Since there thus exists no means of discerning a
        #         # successful from unsuccessful removal attempt, this call is
        #         # currently assumed to always unsuccessfully reduce to a noop.
        #         cell_tri_mesh.remove()
        self._visual.axes.collections = [
            collection
            for collection in self._visual.axes.collections
            if not isinstance(collection, TriMesh)
        ]

        # Return new triangulation meshes for all cells for this time step.
        return self._layer_first_color_mappables()
