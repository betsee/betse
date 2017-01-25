#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector field subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from numpy import ndarray

from betse.lib.numpy import arrays
from betse.science.vector.vectorcls import VectorCells
from betse.util.type.call.memoizers import property_cached
from betse.util.type.types import type_check, NumericTypes, SequenceTypes


# ....................{ CLASSES                            }....................
class VectorField(object):
    '''
    Vector field whose X and Y components are two-dimensional Numpy arrays
    describing a single modelled variable (e.g., total current density)
    spatially situated along one coordinate system (e.g., cell centres, cell
    membrane vertices) for one or more simulation time steps

    This object provides properties efficiently yielding these X and Y
    components in both normalized and nonnormalized forms as well as the
    magnitudes of these components.

    Attributes
    ----------
    _magnitude_factor : NumericTypes
        Factor by which to multiply each magnitude of each vector in this vector
        field, typically to scale vector magnitude to the desired units.
    x : ndarray
        Two-dimensional sequence whose:
        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each X component of each vector in this
          vector field for this time step.
    y : ndarray
        Two-dimensional sequence whose:
        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each Y component of each vector in this
          vector field for this time step.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        x: SequenceTypes,
        y: SequenceTypes,

        # Optional parameters.
        magnitude_factor: NumericTypes = 1,
    ) -> None:
        '''
        Initialize this vector field.

        Parameters
        ----------
        x : ndarray
            Two-dimensional sequence whose:
            * First dimension indexes one or more time steps of this simulation.
            * Second dimension indexes each X component of each vector in this
              vector field for this time step.
        y : ndarray
            Two-dimensional sequence whose:
            * First dimension indexes one or more time steps of this simulation.
            * Second dimension indexes each Y component of each vector in this
              vector field for this time step.
        magnitude_factor : optional[NumericTypes]
            Factor by which to multiply each magnitude of each vector in this
            vector field, typically to scale magnitude to the desired units.
            Defaults to 1, implying no scaling.
        '''

        # Classify the passed sequences as Numpy arrays for efficiency.
        self._x = arrays.from_sequence(x)
        self._y = arrays.from_sequence(y)

        # Classify all remaining parameters.
        self._magnitude_factor = magnitude_factor

    # ..................{ PROPERTIES ~ abstract              }..................
    @property
    def x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each X component of each vector in this
          vector field for this time step.
        '''

        return self._x


    @property
    def y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each Y component of each vector in this
          vector field for this time step.
        '''

        return self._y

    # ..................{ PROPERTIES                         }..................
    @property_cached
    def magnitudes(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each magnitude of each vector in this vector
          field for this time step. For safety, this magnitude is guaranteed to
          be non-zero, avoiding division-by-zero errors elsewhere in the
          codebase when these magnitudes are subsequently divided by (e.g., to
          normalize the X or Y components of this field).

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        # Array of all vector magnitudes computed from the arrays of all
        # vector X and Y components such that each such magnitude is:
        #
        # * Incremented by a negligible positive value approximately equal
        #   to 0, avoiding inevitable division-by-zero errors elsewhere in
        #   the codebase when these magnitudes are subsequently divided by
        #   (e.g., to normalize the X or Y components).
        # * Multiplied by the previously passed factor, typically to scale
        #   magnitudes to the desired units.
        return (
            1e-15 + self._magnitude_factor*np.sqrt(self._x**2 + self._y**2))

    # ..................{ PROPERTIES ~ unit                  }..................
    @property_cached
    def unit_x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary X component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original X component by that vector's original magnitude.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        return self._x / self.magnitudes


    @property_cached
    def unit_y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary Y component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original Y component by that vector's original magnitude.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        return self._y / self.magnitudes

# ....................{ CLASSES ~ cells                    }....................
class VectorFieldCells(object):
    '''
    Cache of various related vector fields whose X and Y components are
    two-dimensional Numpy arrays describing the same underlying data spatially
    situated along different coordinate systems (e.g., cell centres, cell
    membrane vertices) for one or more simulation time steps.

    The input data encapsulated by this cache may be spatially situated at
    either:

    * The centre of each cell in the simulated cluster.
    * The vertex of each cell membrane in the simulated cluster.
    * The midpoint of each cell membrane in the simulated cluster.
    * The centre of each square grid space (in either dimension).

    Each property provided by this vector field (e.g.,
    :meth:`times_cell_centres`) then efficiently interpolates this input data
    from its original coordinate system into the corresponding output data in
    another coordinate system.

    Attributes
    ----------
    _field_args : tuple
        Tuple of all positional arguments to be passed as is to the
        :meth:`VectorCellsABC.__init__` method of each
        :meth:`VectorCellsABC` instance provided by each property of this cache.
    _field_kwargs : dict
        Dictionary of all keyword arguments to be passed as is to the
        :meth:`VectorCellsABC.__init__` method of each
        :meth:`VectorCellsABC` instance provided by each property of this cache.
    _field_x : VectorCells
        Two-dimensional array of the X components of all vectors in this field
        for one or more simulation time steps.
    _field_y : VectorCells
        Two-dimensional array of the Y components of all vectors in this field
        for one or more simulation time steps.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, x: VectorCells, y: VectorCells, *args, **kwargs) -> None:
        '''
        Initialize this vector field cache.

        Parameters
        ----------
        x : VectorCells
            Two-dimensional array of the X components of all vectors in this
            field for one or more simulation time steps.
        y : VectorCells
            Two-dimensional array of the Y components of all vectors in this
            field for one or more simulation time steps.

        All remaining parameters are passed as is to the
        :meth:`VectorCellsABC.__init__` method of each
        :meth:`VectorCellsABC` instance provided by each property of this cache.
        '''

        # Classify all passed parameters.
        self._field_x = x
        self._field_y = y
        self._field_args = args
        self._field_kwargs = kwargs

    # ..................{ PROPERTIES                         }..................
    # Read-only properties, preventing callers from setting these attributes.

    @property_cached
    def times_cells_centre(self) -> VectorField:
        '''
        Vector field whose X and Y components are two-dimensional Numpy arrays
        of data spatially situated at cell centres for one or more time steps.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.

        See Also
        ----------
        :meth:`VectorCells.times_cells_centre`
            Details on the structure of these arrays.
        '''

        return VectorField(
            *self._field_args,
            x=self._field_x.times_cells_centre,
            y=self._field_y.times_cells_centre,
            **self._field_kwargs
        )


    @property_cached
    def times_grids_centre(self) -> VectorField:
        '''
        Vector field whose X and Y components are two-dimensional Numpy arrays
        of data spatially situated at grid space centres for one or more
        time steps.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.

        See Also
        ----------
        :meth:`VectorCells.times_grids_centre`
            Details on the structure of these arrays.
        '''

        return VectorField(
            *self._field_args,
            x=self._field_x.times_grids_centre,
            y=self._field_y.times_grids_centre,
            **self._field_kwargs
        )


    @property_cached
    def times_membranes_midpoint(self) -> VectorField:
        '''
        Vector field whose X and Y components are two-dimensional Numpy arrays
        of data spatially situated at cell membrane midpoints for one or more
        time steps.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.

        See Also
        ----------
        :meth:`VectorCells.times_membranes_midpoint`
            Details on the structure of these arrays.
        '''

        return VectorField(
            *self._field_args,
            x=self._field_x.times_membranes_midpoint,
            y=self._field_y.times_membranes_midpoint,
            **self._field_kwargs
        )


    @property_cached
    def times_membranes_vertex(self) -> VectorField:
        '''
        Vector field whose X and Y components are two-dimensional Numpy arrays
        of data spatially situated at cell membrane vertices for one or more
        time steps.

        For space and time efficiency, the definition of this field is lazily
        deferred to the first read of this property.

        See Also
        ----------
        :meth:`VectorCells.times_membranes_vertex`
            Details on the structure of these arrays.
        '''

        return VectorField(
            *self._field_args,
            x=self._field_x.times_membranes_vertex,
            y=self._field_y.times_membranes_vertex,
            **self._field_kwargs
        )
