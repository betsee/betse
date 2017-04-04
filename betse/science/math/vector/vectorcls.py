#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseVectorException
from betse.science.simulate.simphase import SimPhaseABC
from betse.lib.numpy import arrays
from betse.util.type.call.memoizers import property_cached
from betse.util.type.types import type_check, SequenceOrNoneTypes
from numpy import ndarray
from scipy import interpolate

# ....................{ SUPERCLASSES                       }....................
class VectorCells(object):
    '''
    Cache of various related two-dimensional Numpy arrays describing the same
    underlying data spatially situated along different coordinate systems (e.g.,
    cell centres, cell membrane vertices) for one or more simulation time steps.

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
    _phase : SimPhaseABC
        Current simulation phase.
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
        phase: SimPhaseABC,
        times_cells_centre: SequenceOrNoneTypes = None,
        times_grids_centre: SequenceOrNoneTypes = None,
        times_membranes_midpoint: SequenceOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this vector cache.

        Parameters
        ----------
        phase : SimPhaseABC
            Current simulation phase.
        times_cells_centre : optional[SequenceTypes]
            Two-dimensional sequence of all cell data for a single cell
            membrane-specific modelled variable (e.g., cell electric field
            magnitude) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell, such that each element is
              arbitrary cell data spatially situated at the centre of this cell
              for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_grids_centre`` and ``times_membranes_midpoint`` parameters
            must be non-``None``.
        times_grids_centre : optional[SequenceTypes]
            Two-dimensional sequence of all grid data for a single
            intra- and/or extracellular modelled variable (e.g., total current
            density) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each grid space (in either dimension),
              such that each element is arbitrary grid data spatially
              situated at the centre of this grid space for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_cells_centre`` and ``times_membranes_midpoint`` parameters
            must be non-``None``.
        times_membranes_midpoint : optional[SequenceTypes]
            Two-dimensional sequence of all cell membrane data for a single
            cell membrane-specific modelled variable (e.g., cell membrane
            voltage) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell membrane, such that each
              element is arbitrary cell membrane data spatially situated at the
              midpoint of this membrane for this time step.
            Defaults to ``None``, in which case at least one of the
            ``times_cells_centre`` and ``times_grids_centre`` parameters must be
            non-``None``.

        Raises
        ----------
        BetseVectorException
            If exactly one of the ``times_cells_centre``,
            ``times_grids_centre``, and ``times_membranes_midpoint`` parameters
            is *not* passed.
        '''

        # Initialize our superclass.
        super().__init__()

        # If no sequence was passed, raise an exception.
        if (
            times_cells_centre is None and
            times_grids_centre is None and
            times_membranes_midpoint is None
        ):
            raise BetseVectorException(
                'Parameters "times_cells_centre", "times_grids_centre", and '
                '"times_membranes_midpoint" not passed.')

        # Convert each passed sequence into a Numpy array for efficiency.
        if times_cells_centre is not None:
            times_cells_centre = arrays.from_sequence(times_cells_centre)
        if times_grids_centre is not None:
            times_grids_centre = arrays.from_sequence(times_grids_centre)
        if times_membranes_midpoint is not None:
            times_membranes_midpoint = arrays.from_sequence(
                times_membranes_midpoint)

        # Classify all passed parameters.
        self._phase = phase
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

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell, such that each element is
          arbitrary cell data spatially situated at the centre of this cell for
          this time step.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.
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
            raise BetseVectorException(
                'Properties "times_cells_centre" and '
                '"times_membranes_midpoint" not convertible from '
                'property "times_grids_centre".')

        # Else, remap this vector from cell membrane midpoints to cell centres.
        return self._phase.cells.map_membranes_midpoint_to_cells_centre(
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

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.
        '''

        # If this vector was originally situated at cell centres, return the
        # original array of such data.
        if self._times_membranes_midpoint is not None:
            return self._times_membranes_midpoint

        # Else, remap this vector from cell centres to cell membrane midpoints.
        return self._phase.cells.map_cells_centre_to_membranes_midpoint(
            self.times_cells_centre)


    @property_cached
    def times_membranes_vertex(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell membrane vertex data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane vertex, such that each
          element is arbitrary data spatially situated at this cell membrane
          vertex for this time step.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.
        '''

        return np.dot(
            self.times_membranes_midpoint, self._phase.cells.matrixMap2Verts)


    @property_cached
    def times_regions_centre(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary Voronoi region centre data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each polygonal regions in the Voronoi
          diagram, such that each element is arbitrary data spatially situated
          at the centre of this region for this time step.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.
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

        . First dimension indexes each simulation time step.
        . Second dimension indexes each square grid space (in either dimension),
          such that each element is either:
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

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.
        '''

        # If this vector was originally situated at grid space centres, return
        # the original array of such data.
        if self._times_grids_centre is not None:
            return self._times_grids_centre
        # Else, remap this vector from cell centres to grid space centres.

        # Two-dimensional sequence of input X and Y coordinates to be
        # interpolated from, whose:
        #
        # * First dimension indexes first the X and then Y dimension.
        # * Second dimension indexes each cell, such that each element is the
        #   input coordinate in this dimension of the centre of this cell.
        dimensions_cells_center = (
            self._phase.cells.cell_centres[:, 0],
            self._phase.cells.cell_centres[:, 1])

        # Two-dimensional sequence of output X and Y coordinates to be
        # interpolated into, whose:
        #
        # * First dimension indexes first the X and then Y dimension.
        # * Second dimension indexes each square grid space, such that each
        #   element is the output coordinate in this dimension of the centre of
        #   this grid space.
        dimensions_grids_centre = (self._phase.cells.X, self._phase.cells.Y)

        # Two-dimensional sequence to be returned.
        times_grids_centre = []

        # For each one-dimensional input array of the X or Y components of all
        # vectors situated at the center of each cell for each time step...
        for cells_centre in self.times_cells_centre:
            # One-dimensional output array of the X or Y components of all
            # vectors situated at the centre of each grid space for this time
            # step, interpolated from this input array.
            grids_centre = self._phase.cells.maskECM * interpolate.griddata(
                # Tuple of X and Y coordinates to interpolate from.
                points=dimensions_cells_center,

                # Tuple of X and Y coordinates to interpolate onto.
                xi=dimensions_grids_centre,

                # Array of X or Y input vector components to interpolate.
                values=cells_centre,

                # Machine-readable string specifying the type of interpolation
                # to perform, established by the current configuration.
                method=self._phase.p.interp_type,

                # Default value for all environmental grid points residing
                # outside the convex hull of the cell centers being interpolated
                # from, which are thus non-interpolatable. To nullify all
                # extracellular vectors, this is zero. By default, this is NaN.
                fill_value=0,
            )

            # Append this output array to this output sequence.
            times_grids_centre.append(grids_centre)

        # Return this output sequence, converted from a list into a Numpy array.
        return arrays.from_sequence(times_grids_centre)
