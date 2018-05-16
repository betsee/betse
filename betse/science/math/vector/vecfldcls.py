#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector field subclasses.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.lib.numpy import nparray
from betse.science.math.vector.veccls import VectorCellsCache
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray

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
        '''

        # Classify the passed sequences as Numpy arrays for efficiency.
        self._x = nparray.from_iterable(x)
        self._y = nparray.from_iterable(y)

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

    # ..................{ PROPERTIES ~ magnitudes            }..................
    @property_cached
    def magnitudes(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each (possibly zero) magnitude of each vector
          in this vector field for this time step.

        This array is created only on the first access of this property.

        Caveats
        ----------
        This array may contain zero values and is thus *not* safe for
        general-purpose use as the divisor of division (e.g., to normalize the X
        and Y components of this field). See :meth:`magnitudes_nonzero` for a
        comparable array guaranteed *not* to contain zero values,
        '''

        #FIXME: For efficiency, refactor this to simply call a builtin numpy
        #norm function or method. (Who knew?)

        # Array of all vector magnitudes computed from the arrays of all
        # vector X and Y components.
        return np.sqrt(self._x**2 + self._y**2)


    @property_cached
    def magnitudes_nonzero(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each magnitude of each vector in this vector
          field for this time step. For safety, this magnitude is guaranteed to
          be non-zero, avoiding division-by-zero errors elsewhere in the
          codebase when these magnitudes are subsequently divided by (e.g., to
          normalize the X or Y components of this field). If the unmodified
          value of this magnitude in the :attr:`magnitudes` array is zero, this
          magnitude is replaced with 1.0. This magnitude is indistinguishable
          from magnitudes whose unmodified values are also 1.0 and is thus
          safely usable *only* in contexts where vectors with zero magnitudes
          are ignorable (e.g., as the divisior of division).

        This array is created only on the first access of this property.
        '''

        # Array of all vector magnitudes copied from the original.
        magnitudes_nonzero = self.magnitudes.copy()

        # Substitute all zero magnitudes by 1.0.
        magnitudes_nonzero[magnitudes_nonzero == 0.0] = 1.0

        # Return the resulting array.
        return magnitudes_nonzero

    # ..................{ PROPERTIES ~ unit                  }..................
    @property_cached
    def unit_x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary X component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original X component by that vector's original magnitude.

        This array is created only on the first access of this property.
        '''

        return self._x / self.magnitudes_nonzero


    @property_cached
    def unit_y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary Y component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original Y component by that vector's original magnitude.

        This array is created only on the first access of this property.
        '''

        return self._y / self.magnitudes_nonzero

# ....................{ CLASSES ~ cache                    }....................
class VectorFieldCellsCache(object):
    '''
    Cell cluster vector field cache, persisting all vector fields whose X and Y
    components are two-dimensional Numpy arrays describing the same underlying
    data spatially situated along different coordinate systems (e.g., cell
    centres, cell membrane vertices) for one or more simulation time steps.

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
        :meth:`VectorField.__init__` method of each :class:`VectorField`
        instance returned by each property of this cache.
    _field_kwargs : dict
        Dictionary of all keyword arguments to be passed as is to the
        :meth:`VectorField.__init__` method of each :class:`VectorField`
        instance returned by each property of this cache.
    _field_x : VectorCellsCache
        Two-dimensional array of the X components of all vectors in this field
        for one or more simulation time steps.
    _field_y : VectorCellsCache
        Two-dimensional array of the Y components of all vectors in this field
        for one or more simulation time steps.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        x: VectorCellsCache,
        y: VectorCellsCache,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this vector field cache.

        Parameters
        ----------
        x : VectorCellsCache
            Two-dimensional array of the X components of all vectors in this
            field for one or more simulation time steps.
        y : VectorCellsCache
            Two-dimensional array of the Y components of all vectors in this
            field for one or more simulation time steps.

        All remaining parameters are passed as is to the
        :meth:`VectorField.__init__` method of each :class:`VectorField`
        instance returned by each property of this cache.
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

        This field is created only on the first access of this property.

        See Also
        ----------
        :meth:`VectorCellsCache.times_cells_centre`
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

        This field is created only on the first access of this property.

        See Also
        ----------
        :meth:`VectorCellsCache.times_grids_centre`
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

        This field is created only on the first access of this property.

        See Also
        ----------
        :meth:`VectorCellsCache.times_membranes_midpoint`
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

        This field is created only on the first access of this property.

        See Also
        ----------
        :meth:`VectorCellsCache.times_membranes_vertex`
            Details on the structure of these arrays.
        '''

        return VectorField(
            *self._field_args,
            x=self._field_x.times_membranes_vertex,
            y=self._field_y.times_membranes_vertex,
            **self._field_kwargs
        )
