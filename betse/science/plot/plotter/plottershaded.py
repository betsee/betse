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
    #FIXME: Implement us up.
    @type_check
    def plot(
        self, plot: 'betse.science.plot.PlotCellsABC') -> None:
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

        # Array of arbitrary cell data mapped from the midpoints of each cell
        # membrane onto the vertices of each cell membrane.
        cell_vertex_data = np.dot(plot.cell_data, plot.cells.matrixMap2Verts)

        # If the cell cluster has yet to be triangulated, this *MUST* be the
        # first call to this method. In this case, create this triangulation.
        if self._cell_tri_meshes is None:
            # Array of cell patches at vertices.
            cell_vertex_faces = upscale_data(plot.cells.cell_verts)
        # Else, the cell cluster has already been triangulated. In this case,
        # simply reshade the triangulated mesh for each cell.
        else:
            pass
