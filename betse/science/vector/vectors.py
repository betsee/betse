#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.util.py import references
from betse.lib.numpy import arrays
from betse.util.type.callables import property_cached
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray

# ....................{ SUPERCLASSES                       }....................
class VectorCells(object):
    '''
    Two-dimensional vector of arbitrary cell data for one or more simulation
    time steps, providing properties permitting input data spatially situated
    along one grid system (e.g., cell centres, cell membrane midpoints) to be
    efficiently interpolated along other grid systems.

    The data encapsulated by this vector may be spatially situated at either:

    * The centre of each cell in the simulated cluster.
    * The vertex of each cell membrane in the simulated cluster.
    * The midpoint of each cell membrane in the simulated cluster.

    Attributes
    ----------
    _cells : Cells
        Current cell cluster.
    _times_membranes_midpoint : ndarray
        Two-dimensional Numpy array of all arbitrary cell data for one or more
        simulation time steps spatially situated at cell membrane midpoints,
        returned by the :meth:`times_membranes_midpoint` property.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, times_membranes_midpoint: SequenceTypes) -> None:
        '''
        Initialize this vector.

        Note that this vector is _not_ safely usable until after the
        :meth:`prep` method has been externally called.

        Parameters
        ----------
        times_membranes_midpoint : SequenceTypes
            Two-dimensional sequence of all cell membrane data for a single
            cell membrane-specific modelled variable (e.g., cell membrane
            voltage) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell membrane, such that each
              element is arbitrary cell membrane data spatially situated at the
              midpoint of this membrane for this time step.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify the passed sequence as a Numpy array for efficiency.
        self._times_membranes_midpoint = arrays.from_sequence(
            times_membranes_midpoint)


    @type_check
    def prep(self, cells: 'betse.science.cells.Cells') -> None:
        '''
        Prepare this vector for subsequent use.

        Note that this vector is _not_ safely usable until after this method has
        been externally called.

        Parameters
        ----------
        cells : Cells
            Current cell cluster.
        '''

        # Classify this cell cluster with a weak rather than strong (the
        # default) reference, thereby avoiding circular references and the
        # resulting complications thereof (e.g., increased memory overhead).
        # Since the parent object necessarily lives significantly longer than
        # this vector, no complications arise. Ergo, this attribute *ALWAYS*
        # yields the expected object (rather than non-deterministically yielding
        # None if this object is unexpectedly garbage-collected).
        self._cells = references.proxy_weak(cells)

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from setting these attributes.

    @property_cached
    def times_cells_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell, such that each element is
          arbitrary cell data spatially situated at the the centre of this cell
          for this time step.
        '''

        return self._cells.map_membranes_midpoint_to_cells_centre(
            self._times_membranes_midpoint)


    @property
    def times_membranes_midpoint(self) -> ndarray:
        '''
        Two-dimensional sequence of all arbitrary cell membrane data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane, such that each element is
          arbitrary cell membrane data spatially situated at the midpoint of
          this membrane for this time step.
        '''

        return self._times_membranes_midpoint


    #FIXME: Refactor the following methods to internally defer to the
    #corresponding Cells.map_*() methods.

    @property_cached
    def times_membranes_vertex(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell membrane vertex data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane vertex, such that each
          element is arbitrary data spatially situated at this cell membrane
          vertex for this time step.
        '''

        #FIXME: Consider generalizing this logic into a new public "Cells"
        #method ala the Cells.map_mems_midpoint_to_cells_centre()
        #method called above.

        return np.dot(
            self._times_membranes_midpoint, self._cells.matrixMap2Verts)


    @property_cached
    def times_regions_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary Voronoi region centre data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each polygonal regions in the Voronoi
          diagram, such that each element is arbitrary data spatially situated
          at the centre of this region for this time step.
        '''

        #FIXME: Consider generalizing this logic into a new public "Cells"
        #method ala the Cells.map_mems_midpoint_to_cells_centre()
        #method called above.

        # Array to be returned, zeroed to the desired shape.
        times_regions_centre = np.zeros((
            len(self.times_cells_centre),
            len(self._cells.voronoi_centres)))

        # Map cell- to region-centred data for all time steps.
        times_regions_centre[
            :, self._cells.cell_to_grid] = self.times_cells_centre

        # Return this array for subsequent caching.
        return times_regions_centre
