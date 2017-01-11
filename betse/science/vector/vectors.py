#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseVectorException
from betse.lib.numpy import arrays
from betse.util.py import references
from betse.util.type.callables import property_cached
from betse.util.type.types import type_check, SequenceOrNoneTypes
from numpy import ndarray

# ....................{ SUPERCLASSES                       }....................
class VectorCells(object):
    '''
    Two-dimensional vector of arbitrary cell data for one or more simulation
    time steps.

    The input data encapsulated by this vector may be spatially situated at
    either:

    * The centre of each cell in the simulated cluster.
    * The vertex of each cell membrane in the simulated cluster.
    * The midpoint of each cell membrane in the simulated cluster.
    * The centre of each square grid space (in either dimension).

    Each property provided by this vector (e.g., :meth:`times_cell_centres`)
    then efficiently interpolates this input data from its original coordinate
    system into the corresponding output data in another coordinate system.

    Attributes
    ----------
    _cells : Cells
        Current cell cluster.
    _p : Parameters
        Current simulation configuration.
    _times_cell_centres : ndarray
        Two-dimensional Numpy array of all arbitrary cell data for one or more
        simulation time steps spatially situated at cell centres, returned by
        the :meth:`times_cells_centre` property.
    _times_membranes_midpoint : ndarray
        Two-dimensional Numpy array of all arbitrary cell membrane data for one
        or more simulation time steps spatially situated at cell membrane
        midpoints, returned by the :meth:`times_membranes_midpoint` property.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        times_cells_centre: SequenceOrNoneTypes = None,
        times_membranes_midpoint: SequenceOrNoneTypes = None,
        cells: 'betse.science.cells.Cells' = None,

        #FIXME: Uncomment this *AFTER* improving the @type_check decorator to
        #handle strings embedded in annotation tuples.
        # cells: ('betse.science.cells.Cells', NoneType) = None,
    ) -> None:
        '''
        Initialize this vector.

        Note that this vector is _not_ safely usable until after the
        :meth:`prep` method has been externally called.

        Parameters
        ----------
        p : Parameters
            Current simulation configuration.
        cells : betse.science.cells.Cells
            Current cell cluster.
        times_cells_centre : optional[SequenceTypes]
            Two-dimensional sequence of all cell data for a single cell
            membrane-specific modelled variable (e.g., cell electric field
            magnitude) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell, such that each element is
              arbitrary cell data spatially situated at the centre of this cell
              for this time step.
            The `times_membranes_midpoint` parameter _must_ be:
            * `None` if this parameter is non-`None`.
            * Non-`None` if this parameter is `None`.
            In other words, exactly one of this and the
            `times_membranes_midpoint` parameters _must_ be passed. Defaults to
            `None`.
        times_membranes_midpoint : optional[SequenceTypes]
            Two-dimensional sequence of all cell membrane data for a single
            cell membrane-specific modelled variable (e.g., cell membrane
            voltage) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell membrane, such that each
              element is arbitrary cell membrane data spatially situated at the
              midpoint of this membrane for this time step.
            The `times_cells_centre` parameter _must_ be:
            * `None` if this parameter is non-`None`.
            * Non-`None` if this parameter is `None`.
            In other words, exactly one of this and the `times_cells_centre`
            parameters _must_ be passed. Defaults to `None`.

        Raises
        ----------
        BetseVectorException
            If exactly one of the `times_cells_centre` and
            `times_membranes_midpoint` parameters is *not* passed.
        '''

        # Initialize our superclass.
        super().__init__()

        # If no sequence was passed, raise an exception.
        if times_cells_centre is None and times_membranes_midpoint is None:
            raise BetseVectorException(
                'Parameters "times_cells_centre" and '
                '"times_membranes_midpoint" not passed.')

        # If both sequences were passed, raise an exception.
        if (
            times_cells_centre is not None and
            times_membranes_midpoint is not None
        ):
            raise BetseVectorException(
                'Parameters "times_cells_centre" and '
                '"times_membranes_midpoint" both passed.')

        # Default all instance variables.
        self._cells = None
        self._times_cells_centre = None
        self._times_membranes_midpoint = None

        # Classify each passed sequence as a Numpy array for efficiency.
        if times_cells_centre is not None:
            self._times_cells_centre = arrays.from_sequence(times_cells_centre)
        if times_membranes_midpoint is not None:
            self._times_membranes_midpoint = arrays.from_sequence(
                times_membranes_midpoint)

        # If passed a cell cluster, prepare this vector with this cluster.
        if cells is not None:
            self.prep(cells=cells)


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

    # ..................{ PROPERTIES                         }..................
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

        # If this vector was originally situated at cell centres, return the
        # original array of such data.
        if self._times_cells_centre is not None:
            return self._times_cells_centre

        # Else, remap this vector from cell membrane midpoints to cell centres.
        return self._cells.map_membranes_midpoint_to_cells_centre(
            self.times_membranes_midpoint)


    @property_cached
    def times_membranes_midpoint(self) -> ndarray:
        '''
        Two-dimensional sequence of all arbitrary cell membrane data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane, such that each element is
          arbitrary cell membrane data spatially situated at the midpoint of
          this membrane for this time step.
        '''

        # If this vector was originally situated at cell centres, return the
        # original array of such data.
        if self._times_membranes_midpoint is not None:
            return self._times_membranes_midpoint

        # Else, remap this vector from cell centres to cell membrane midpoints.
        return self._cells.map_cells_centre_to_membranes_midpoint(
            self.times_cells_centre)


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
            self.times_membranes_midpoint, self._cells.matrixMap2Verts)


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

    # ..................{ PROPERTIES ~ grid                  }..................
    # Read-only properties, preventing callers from setting these attributes.

    #FIXME: Generalize from the current grid interpolation implemented by the
    #VectorFieldCurrentIntra._get_component() method. Note that, to do so, we'll
    #need to require that an additional "p" parameter be passed to the
    #__init__() method. For sanity, both this and the existing "cells" parameter
    #should be made mandatory. To avoid future refactoring pain, ensure we do
    #this *NOW* before things spiral out-of-hand elsewhere in the codebase.

    @property_cached
    def times_grid_centres(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary grid data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each square grid space (in either dimension),
          such that each element is arbitrary grid data spatially situated at
          the the centre of this grid space for this time step.
        '''

        pass
