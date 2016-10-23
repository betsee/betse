#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying streamlines onto the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
import weakref
from betse.exceptions import BetseSimConfigException
from betse.science.vector.fieldabc import VectorFieldABC
from betse.util.type.types import type_check, SequenceTypes
from scipy import interpolate

# ....................{ SUPERCLASSES                       }....................
class VectorFieldCurrentABC(VectorFieldABC):
    '''
    Abstract base class of all current density vector field subclasses.

    Each subclass of this class is a vector field of current densities of all
    intracellular and/or extracellular spaces for all time steps of the current
    simulation.

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
        sim: 'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p: 'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Initialize this vector field.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        '''

        # Initialize the superclass.
        super().__init__(
            # Upscale the magnitude of each current density vector, yielding
            # magnitude in units of uA/cm^2.
            magnitude_factor=100,
        )

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # objects necessarily live significantly longer than this plot, no
        # complications arise. Ergo, these attributes *ALWAYS* yield these
        # objects rather than non-deterministically yielding "None" if these
        # objects are unexpectedly garbage-collected.
        self._sim = weakref.proxy(sim)
        self._cells = weakref.proxy(cells)
        self._p = weakref.proxy(p)

# ....................{ SUBCLASSES                         }....................
class VectorFieldCurrentIntraExtra(VectorFieldCurrentABC):
    '''
    Vector field of the current density of all intracellular and extracellular
    spaces for all time steps of the current simulation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize this vector field, passing all passed arguments as is to the
        superclass initializer.

        Raises
        ----------
        BetseSimConfigException
            If extracellular spaces are disabled.
        '''

        # Initialize the superclass with all passed arguments.
        super().__init__(*args, **kwargs)

        # If extracellular spaces are disabled, raise an exception.
        if not self._p.sim_ECM:
            raise BetseSimConfigException('Extracellular spaces disabled.')

    # ..................{ SUPERCLASS                         }..................
    @property
    def x(self) -> None:
        return self._sim.I_tot_x_time


    @property
    def y(self) -> None:
        return self._sim.I_tot_y_time


class VectorFieldCurrentIntra(VectorFieldCurrentABC):
    '''
    Vector field of the current density of all intracellular spaces (e.g., gap
    junctions) for all time steps of the current simulation.

    Attributes
    ----------
    _x : ndarray
        Two-dimensional Numpy array as documented by the :meth:`x` property if
        that property has been read at least once for this instance _or_ `None`
        otherwise (in which case this attribute is defined on the first read of
        that property).
    _y : ndarray
        Two-dimensional Numpy array as documented by the :meth:`y` property if
        that property has been read at least once for this instance _or_ `None`
        otherwise (in which case this attribute is defined on the first read of
        that property).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize this vector field, passing all passed arguments as is to the
        superclass initializer.
        '''

        # Initialize the superclass with all passed arguments.
        super().__init__(*args, **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property
    def x(self) -> None:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each X component of each intracellular
          current density vector in this vector field for this time step.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        # If this property has yet to be read, define and cache this array.
        if self._x is None:
            self._init_components()

        # Return the previously cached array.
        return self._x


    @property
    def y(self) -> None:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes one or more time steps of this simulation.
        * Second dimension indexes each Y component of each intracellular
          current density vector in this vector field for this time step.

        For space and time efficiency, the definition of this array is lazily
        deferred to the first read of this property.
        '''

        # If this property has yet to be read, define and cache this array.
        if self._y is None:
            self._init_components()

        # Return the previously cached array.
        return self._y

    # ..................{ SUPERCLASS                         }..................
    def _init_components(self) -> None:
        '''
        Interpolate the X and Y components of all intracellular current density
        vectors in this vector field for all time steps from the center of each
        cell to the X/Y grid of this cell cluster's environment.

        In the case of the :class:`VectorFieldCurrentIntraExtra` subclass, the
        current density of all intracellular and extracellular spaces is
        precomputed by the :class:`Simulator` class and trivially usable as is.
        In the case of this subclass, however, the current density of only
        intracellular spaces is _not_ precomputed elsewhere and hence must be
        computed here.
        '''

        # Two-dimensional tuple of the X and Y components of all cell centers.
        dimensions_cells_center = (
            self._visual.cells.cell_centres[:, 0],
            self._visual.cells.cell_centres[:, 1])

        # Two-dimensional tuple of the X and Y components of all grid points.
        dimensions_grid_points = (self._visual.cells.X, self._visual.cells.Y)

        # Two-dimensional lists of the X and Y components of all current density
        # vectors over all time steps.
        self._time_currents_x = self._get_component(
            times_currents_center=self._visual.sim.I_cell_x_time,
            dimensions_cells_center=dimensions_cells_center,
            dimensions_grid_points=dimensions_grid_points,
        )
        self._time_currents_y = self._prep_currents_component(
            times_currents_center=self._visual.sim.I_cell_y_time,
            dimensions_cells_center=dimensions_cells_center,
            dimensions_grid_points=dimensions_grid_points,
        )


    @type_check
    def _get_component(
        self,
        times_currents_center:   SequenceTypes,
        dimensions_cells_center: SequenceTypes,
        dimensions_grid_points:  SequenceTypes,
    ) -> SequenceTypes:
        '''
        Interpolate the passed array of the X and Y components of all
        intracellular current density vectors in this vector field for all time
        steps from the passed cell centers to the passed X/Y grid of this cell
        cluster's environment.

        Parameters
        ----------
        times_currents_center : SequenceTypes
            Two-dimensional sequence whose:
            * First dimension indexes each time step of this simulation.
            * Second dimension indexes the X or Y components of all
              intracellular current density vectors computed at the center of
              each cell for the corresponding time step.
        dimensions_cells_center : SequenceTypes
            Two-dimensional sequence, whose:
            * First dimension indexes first X and then Y dimensions.
            * Second dimension indexes cells, whose elements are the coordinates
              in that X or Y dimension of the centers of those cells.
        dimensions_grid_points : SequenceTypes
            Two-dimensional sequence, whose:
            * First dimension indexes first X and then Y dimensions.
            * Second dimension indexes the coordinates in that X or Y dimension
              of each grid point.

        Returns
        ----------
        SequenceTypes
            Two-dimensional sequence whose:
            * First dimension indexes each time step of this simulation.
            * Second dimension indexes the X or Y components of all current
              density vectors computed for the corresponding time step.
        '''

        # Output two-dimensional sequence as documented above.
        time_currents_grid = []

        # For each one-dimensional array of the X or Y components of all
        # intracellular current density vectors computed at the center of each
        # cell for each time step...
        for currents_center in times_currents_center:
            # One-dimensional array of the X or Y components of all
            # intracellular current density vectors defined on the X/Y grid of
            # this cell cluster's environment for this time step, interpolated
            # from the corresponding components defined at cell centers.
            currents_grid = self._cells.maskECM * interpolate.griddata(
                # Tuple of X and Y coordinates to interpolate from.
                points=dimensions_cells_center,

                # Tuple of X and Y coordinates to interpolate onto.
                xi=dimensions_grid_points,

                # Array of X or Y current vector components to interpolate.
                values=currents_center,

                # Machine-readable string specifying the type of interpolation
                # to perform, established by the current configuration.
                method=self._p.interp_type,

                # Default value for all environmental grid points residing
                # outside the convex hull of the cell centers being interpolated
                # from, which are thus non-interpolatable. To nullify all
                # extracellular current density vectors. this is zero. By
                # default, this is NaN.
                fill_value=0,
            )

            # Append this grid-interpolated to this output sequence, providing
            # the X or Y components of all current vectors for this time step.
            time_currents_grid.append(currents_grid)

        # Return this output sequence.
        return time_currents_grid
