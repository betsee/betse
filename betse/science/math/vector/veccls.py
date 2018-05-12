#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseSimVectorException
from betse.lib.numpy import nparray
from betse.science.math.cache.cacheabc import SimPhaseCacheABC
from betse.util.type import sequences
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, IterableOrNoneTypes
from numpy import ndarray

# ....................{ SUPERCLASSES                       }....................
class VectorCellsCache(SimPhaseCacheABC):
    '''
    Cell cluster vector cache, persisting all two-dimensional Numpy arrays
    describing the same underlying data spatially situated along different
    coordinate systems (e.g., cell centres, cell membrane vertices) for one or
    more simulation time steps.

    The input data encapsulated by this cache may be spatially situated at
    either:

    * The centre of each cell in the simulated cluster.
    * The vertex of each cell membrane in the simulated cluster.
    * The midpoint of each cell membrane in the simulated cluster.
    * The centre of each square grid space (in either dimension).

    Each property provided by this cache (e.g., :meth:`times_cell_centres`) then
    efficiently interpolates this input data from its original coordinate system
    into the corresponding output data in another coordinate system.

    Attributes
    ----------
    _times_cell_centres : ndarray
        Two-dimensional Numpy array of all arbitrary cell data for one or more
        simulation time steps spatially situated at cell centres, returned by
        the :meth:`times_cells_centre` property.
    _times_membranes_midpoint : ndarray
        Two-dimensional Numpy array of all arbitrary cell membrane data for one
        or more simulation time steps spatially situated at cell membrane
        midpoints, returned by the :meth:`times_membranes_midpoint` property.
    _times_grids_centre : ndarray
        Two-dimensional Numpy array of all arbitrary grid data for one or more
        simulation time steps spatially situated at grid space centres, returned
        by the :meth:`times_grids_centre` property.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        times_cells_centre: IterableOrNoneTypes = None,
        times_grids_centre: IterableOrNoneTypes = None,
        times_membranes_midpoint: IterableOrNoneTypes = None,
        **kwargs
    ) -> None:
        '''
        Initialize this cache.

        Parameters
        ----------
        times_cells_centre : optional[IterableTypes]
            Two-dimensional iterable of all cell data for a single cell
            membrane-specific modelled variable (e.g., cell electric field
            magnitude) for all simulation time steps, whose:
            . First dimension indexes each sampled time step.
            . Second dimension indexes each cell, such that each element is
              arbitrary cell data spatially situated at the centre of this cell
              for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_grids_centre`` and ``times_membranes_midpoint`` parameters
            must be non-``None``.
        times_grids_centre : optional[IterableTypes]
            Two-dimensional iterable of all grid data for a single
            intra- and/or extracellular modelled variable (e.g., total current
            density) for all simulation time steps, whose:
            . First dimension indexes each sampled time step.
            . Second dimension indexes each grid space (in either dimension),
              such that each element is arbitrary grid data spatially
              situated at the centre of this grid space for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_cells_centre`` and ``times_membranes_midpoint`` parameters
            must be non-``None``.
        times_membranes_midpoint : optional[IterableTypes]
            Two-dimensional iterable of all cell membrane data for a single
            cell membrane-specific modelled variable (e.g., cell membrane
            voltage) for all simulation time steps, whose:
            . First dimension indexes each sampled time step.
            . Second dimension indexes each cell membrane, such that each
              element is arbitrary cell membrane data spatially situated at the
              midpoint of this membrane for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_cells_centre`` and ``times_grids_centre`` parameters must be
            non-``None``.

        All remaining keyword arguments are passed as is to the superclass
        :meth:`SimPhaseCacheABC.__init__` method.

        Raises
        ----------
        BetseSimVectorException
            If exactly one of the ``times_cells_centre``,
            ``times_grids_centre``, and ``times_membranes_midpoint`` parameters
            is *not* passed.
        BetseSequenceException
            If any of the ``times_cells_centre``, ``times_grids_centre``, and
            ``times_membranes_midpoint`` parameters that are passed are empty
            (i.e., sequences of length 0).
        '''

        # Initialize our superclass.
        super().__init__(**kwargs)

        # If no iterable was passed, raise an exception.
        if (
            times_cells_centre is None and
            times_grids_centre is None and
            times_membranes_midpoint is None
        ):
            raise BetseSimVectorException(
                'Parameters "times_cells_centre", "times_grids_centre", and '
                '"times_membranes_midpoint" not passed.')

        # Convert each passed iterable into a Numpy array for efficiency,
        # raising an exception if any such iterable is empty.
        if times_cells_centre is not None:
            times_cells_centre = nparray.from_iterable(times_cells_centre)
            sequences.die_if_empty(
                times_cells_centre,
                label='Parameter "times_cells_centre"')
        if times_grids_centre is not None:
            times_grids_centre = nparray.from_iterable(times_grids_centre)
            sequences.die_if_empty(
                times_grids_centre,
                label='Parameter "times_grids_centre"')
        if times_membranes_midpoint is not None:
            times_membranes_midpoint = nparray.from_iterable(
                times_membranes_midpoint)
            sequences.die_if_empty(
                times_membranes_midpoint,
                label='Parameter "times_membranes_midpoint"')

        # Classify all passed parameters.
        self._times_cells_centre = times_cells_centre
        self._times_grids_centre = times_grids_centre
        self._times_membranes_midpoint = times_membranes_midpoint

    # ..................{ PROPERTIES                         }..................
    # Read-only properties, preventing callers from setting these attributes.

    @property_cached
    def times_cells_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell data for all
        simulation time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell, such that each element is
           arbitrary cell data spatially situated at the centre of this cell for
           this time step.

        This array is created only on the first access of this property.
        '''

        # If this vector was originally situated at cell centres, return the
        # original array of such data.
        if self._times_cells_centre is not None:
            return self._times_cells_centre

        # If this vector is *NOT* situated at cell membrane midpoints (e.g., is
        # situated at grid space centres), raise an exception. In this case,
        # permitting the times_membranes_midpoint() method to be called below
        # would recall this method, which would recall that method, and so on
        # until the stack is exhausted, halting the current process.
        if self._times_membranes_midpoint is None:
            raise BetseSimVectorException(
                'Properties "times_cells_centre" and '
                '"times_membranes_midpoint" not convertible from '
                'property "times_grids_centre".')

        # Else, remap this vector from cell membrane midpoints to cell centres.
        return self._phase.cells.map_membranes_midpoint_to_cells_centre(
            membranes_midpoint_data=self.times_membranes_midpoint)


    @property_cached
    def times_membranes_midpoint(self) -> ndarray:
        '''
        Two-dimensional sequence of all arbitrary cell membrane data for all
        simulation time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell membrane, such that each element is
           arbitrary cell membrane data spatially situated at the midpoint of
           this membrane for this time step.

        This array is created only on the first access of this property.
        '''

        # If this vector was originally situated at cell membrane midpoints,
        # return the original array of such data.
        if self._times_membranes_midpoint is not None:
            return self._times_membranes_midpoint

        # Else, remap this vector from cell centres to cell membrane midpoints.
        return self._phase.cells.map_cells_centre_to_membranes_midpoint(
            cells_centre_data=self.times_cells_centre)


    @property_cached
    def times_membranes_vertex(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell membrane vertex data
        for all simulation time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell membrane vertex, such that each
           element is arbitrary data spatially situated at this cell membrane
           vertex for this time step.

        This array is created only on the first access of this property.
        '''

        return np.dot(
            self.times_membranes_midpoint, self._phase.cells.matrixMap2Verts)


    @property_cached
    def times_regions_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary Voronoi region centre data
        for all simulation time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each polygonal regions in the Voronoi
           diagram, such that each element is arbitrary data spatially situated
           at the centre of this region for this time step.

        This array is created only on the first access of this property.
        '''

        # Array to be returned, zeroed to the desired shape.
        times_regions_centre = np.zeros((
            len(self.times_cells_centre),
            len(self._phase.cells.voronoi_centres)))

        # Map cell- to region-centred data for all time steps.
        times_regions_centre[
            :, self._phase.cells.cell_to_grid] = self.times_cells_centre

        # Return this array for subsequent caching.
        return times_regions_centre


    #FIXME: Is this actually grid point vertices rather than grid space centres?
    #It doesn't particularly matter in terms of implementation (which clearly
    #works), but it would be useful to eliminate incorrectness in terminology.
    @property_cached
    def times_grids_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary gridded cell data for all
        simulation time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each square environmental grid space (in
           either dimension), such that each element is either:
           * If this vector was initialized with the :meth:`times_grids_centre`
             parameter, arbitrary grid data spatially situated at the centre of
             this grid space for this time step. In this case, no interpolation
             is required and the original :meth:`times_grids_centre` array is
             returned as is. This edge case preserves all extracellular data for
             environmental grid spaces.
           * Else:
             * If this grid space resides *inside* the convex hull of this cell
               cluster, arbitrary grid data spatially interpolated from the
               centres of all cells whose membranes overlap this grid space onto
               the centre of this grid space for this time step.
             * Else, 0. In this case, this is an environmental grid space
               residing *outside* the convex hull of this cell cluster. Since
               this vector was initialized with a cell-centric array rather than
               the `times_grids_centre` parameter, this vector contains no
               extracellular data to interpolate environmental grid spaces from.

        This array is created only on the first access of this property.
        '''

        # If this vector was originally situated at grid space centres, return
        # the original array of such data.
        if self._times_grids_centre is not None:
            # print('Reusing cached "times_grids_centre"...')
            return self._times_grids_centre
        # print('Caching new "times_grids_centre"...')

        # Else, remap this vector from cell centres to grid space centres.
        return self._phase.cells.map_cells_centre_to_grids_centre(
            cells_centre_data=self.times_cells_centre,
            interp_method=self._phase.p.interp_type,
        )
