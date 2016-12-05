#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all vector field subclasses.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty

import numpy as np
from betse.util.py import references
from betse.util.type.obj.objs import property_cached
from betse.util.type.types import type_check, NumericTypes
from numpy import ndarray


# ....................{ CLASSES                            }....................
class VectorFieldABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all vector field subclasses.

    Each subclass of this class provides properties efficiently yielding all
    X and Y components and magnitudes of the vector field of a single modelled
    variable (e.g., current density) for one or more simulation time steps.

    Attributes
    ----------
    _magnitude_factor : NumericTypes
        Factor by which to multiply each magnitude of each vector in this vector
        field, typically to scale magnitude to the desired units.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, magnitude_factor: NumericTypes = 1) -> None:
        '''
        Initialize this vector field.

        Parameters
        ----------
        magnitude_factor : NumericTypes
            Factor by which to multiply each magnitude of each vector in this
            vector field, typically to scale magnitude to the desired units.
            Defaults to 1, implying no scaling.
        '''

        # Classify the passed parameters.
        self._magnitude_factor = magnitude_factor

    # ..................{ PROPERTIES ~ abstract              }..................
    @abstractproperty
    def x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each X component of each vector in this
          vector field for this time step.
        '''

        pass


    @abstractproperty
    def y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each Y component of each vector in this
          vector field for this time step.
        '''

        pass

    # ..................{ PROPERTIES ~ concrete              }..................
    @property_cached
    def magnitudes(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each magnitude of each vector in this vector
          field for this time step.

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
            1e-15 + self._magnitude_factor*np.sqrt(self.x**2 + self.y**2))


    @property_cached
    def x_unit(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary X component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original X component by that vector's original magnitude.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        return self.x / self._magnitudes


    @property_cached
    def y_unit(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each unitary Y component of each unit vector
          in this vector field for this time step, produced by dividing that
          vector's original Y component by that vector's original magnitude.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        return self.y / self._magnitudes


class VectorFieldSimulatedABC(VectorFieldABC):
    '''
    Abstract base class of all **simulated vector field subclasses** (i.e.,
    vector fields simulated by the current simulation and hence attributes of
    the current :class:`Simulator` instance).

    Attributes
    ----------
    _cells : Cells
        Current cell cluster.
    _p : Parameters
        Current simulation configuration.
    _sim : Simulator
        Current simulation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        *args, **kwargs
    ) -> None:
        '''
        Initialize this simulated vector field.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.

        All remaining parameters are passed as is to the superclass.
        '''

        # Initialize our superclass with all remaining parameters.
        super().__init__(*args, **kwargs)

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # objects necessarily live significantly longer than this plot, no
        # complications arise. Ergo, these attributes *ALWAYS* yield these
        # objects rather than non-deterministically yielding "None" if these
        # objects are unexpectedly garbage-collected.
        self._sim = references.proxy_weak(sim)
        self._cells = references.proxy_weak(cells)
        self._p = references.proxy_weak(p)
