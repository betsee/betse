#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying streamlines onto the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimConfigException
from betse.lib.numpy import arrays
from betse.science.vector.fieldabc import VectorFieldABC
from betse.util.py import references
from betse.util.type.objects import property_cached
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray
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
    def __init__(
        self,
        sim: 'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p: 'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Initialize this current density vector field.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        '''

        # Initialize our superclass.
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
        self._sim = references.proxy_weak(sim)
        self._cells = references.proxy_weak(cells)
        self._p = references.proxy_weak(p)

# ....................{ SUBCLASSES                         }....................
class VectorFieldCurrentIntraExtra(VectorFieldCurrentABC):
    '''
    Vector field of the current densities of all intracellular and extracellular
    spaces for all time steps of the current simulation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # If extracellular spaces are disabled, raise an exception.
        if not self._p.sim_ECM:
            raise BetseSimConfigException('Extracellular spaces disabled.')

    # ..................{ SUBCLASS                           }..................
    # Reuse the X and Y components already computed for this simulation,
    # converted from lists to Numpy arrays.

    @property_cached
    def x(self) -> ndarray:
        return arrays.from_sequence(self._sim.I_tot_x_time)


    @property_cached
    def y(self) -> ndarray:
        return arrays.from_sequence(self._sim.I_tot_y_time)


class VectorFieldCurrentIntra(VectorFieldCurrentABC):
    '''
    Vector field of the current densities of all intracellular spaces (e.g., gap
    junctions) for all time steps of the current simulation.
    '''

    # ..................{ SUPERCLASS                         }..................
    @property_cached
    def x(self) -> ndarray:
        return self._get_component(self._sim.I_cell_x_time)


    @property_cached
    def y(self) -> ndarray:
        return self._get_component(self._sim.I_cell_y_time)

    # ..................{ GETTERS                            }..................
    @type_check
    def _get_component(self, times_currents_centered: SequenceTypes) -> ndarray:
        '''
        Interpolate the passed array of the X and Y components of all
        intracellular current density vectors in this vector field for all time
        steps from the passed cell centers to the passed X/Y grid of this cell
        cluster's environment.

        In the case of the :class:`VectorFieldCurrentIntraExtra` subclass, the
        current density of all intracellular and extracellular spaces is
        precomputed by the :class:`Simulator` class and trivially usable as is.
        In the case of this subclass, however, the current density of only
        intracellular spaces is _not_ precomputed elsewhere and hence must be
        computed here.

        Parameters
        ----------
        times_currents_center : SequenceTypes
            Two-dimensional sequence whose:
            * First dimension indexes each time step of this simulation.
            * Second dimension indexes the X or Y components of all
              intracellular current density vectors computed at the center of
              each cell for the corresponding time step.

        Returns
        ----------
        ndarray
            Two-dimensional Numpy array whose:
            * First dimension indexes each time step of this simulation.
            * Second dimension indexes the X or Y components of all current
              density vectors computed for the corresponding time step.
        '''

        # Output two-dimensional sequence as documented above.
        time_currents_gridded = []

        # Tuple of X and Y coordinates to interpolate from and onto
        # (respectively), localized for premature optimization.
        dimensions_cells_center = self._dimensions_cells_center

        # Tuple of X and Y coordinates to interpolate onto.
        dimensions_grid_points = self._dimensions_grid_points

        # For each one-dimensional array of the X or Y components of all
        # intracellular current density vectors computed at the center of each
        # cell for each time step...
        for currents_centered in times_currents_centered:
            # One-dimensional array of the X or Y components of all
            # intracellular current density vectors defined on the X/Y grid of
            # this cell cluster's environment for this time step, interpolated
            # from the corresponding components defined at cell centers.
            currents_gridded = self._cells.maskECM * interpolate.griddata(
                # Tuple of X and Y coordinates to interpolate from.
                points=dimensions_cells_center,

                # Tuple of X and Y coordinates to interpolate onto.
                xi=dimensions_grid_points,

                # Array of X or Y current vector components to interpolate.
                values=currents_centered,

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
            time_currents_gridded.append(currents_gridded)

        # Return this output sequence, converted from a list to Numpy array.
        return arrays.from_sequence(time_currents_gridded)

    # ..................{ PROPERTIES                         }..................
    @property_cached
    def _dimensions_cells_center(self) -> tuple:
        '''
        Two-dimensional tuple of the X and Y components of all cell centers.

        Two-dimensional sequence, whose:
        * First dimension indexes first X and then Y dimensions.
        * Second dimension indexes cells, whose elements are the coordinates
          in that X or Y dimension of the centers of those cells.
        '''

        return (self._cells.cell_centres[:, 0], self._cells.cell_centres[:, 1])


    @property_cached
    def _dimensions_grid_points(self) -> tuple:
        '''
        Two-dimensional sequence, whose:

        * First dimension indexes first X and then Y dimensions.
        * Second dimension indexes the coordinates in that X or Y dimension of
          each uniformally spaced grid point for both this cell cluster and its
          external environment.
        '''

        return (self._cells.X, self._cells.Y)
