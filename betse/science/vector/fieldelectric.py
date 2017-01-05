#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Vector field classes describing electromagnetic phenomena, including both
intracellular and extracellular simulated fields and current density fields.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.numpy import arrays
from betse.science.vector.fieldabc import VectorFieldSimmedABC
from betse.util.type.callables import property_cached
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray
from scipy import interpolate

# ....................{ CONSTANTS                          }....................
_MAGNITUDE_FACTOR_CURRENT = 100
'''
Factor by which to multiply each magnitude of each vector in each current
density vector field, yielding magnitude in units of uA/cm^2.
'''

# ....................{ ELECTRIC FIELD                     }....................

# ....................{ CURRENT DENSITY ~ intra-extra      }....................
class VectorFieldCurrentIntraExtra(VectorFieldSimmedABC):
    '''
    Vector field of the current densities of all intracellular and extracellular
    spaces spatially situated at grid space centres for all time steps of the
    current simulation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass to...
        super().__init__(
            # Require simulation of extracellular spaces.
            is_ecm_needed=True,

            # Upscale vector magnitudes by this multiplicative factor.
            magnitude_factor=_MAGNITUDE_FACTOR_CURRENT,

            # Pass all remaining parameters as is.
            *args, **kwargs)

    # ..................{ SUBCLASS                           }..................
    # Reuse the X and Y components already computed for this simulation,
    # converted from lists to Numpy arrays.

    @property_cached
    def x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the X
          component of the **total current density vector** (i.e., vector of
          both intra- _and_ extracellular current densities) spatially situated
          at the center of each grid space for this time step.
        '''

        return arrays.from_sequence(self._sim.I_tot_x_time)


    @property_cached
    def y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the Y
          component of the **total current density vector** (i.e., vector of
          both intra- _and_ extracellular current densities) spatially situated
          at the center of each grid space for this time step.
        '''

        return arrays.from_sequence(self._sim.I_tot_y_time)

# ....................{ CURRENT DENSITY ~ intra            }....................
class VectorFieldCurrentIntra(VectorFieldSimmedABC):
    '''
    Vector field of the current densities of all intracellular spaces (e.g., gap
    junctions) spatially situated at grid space centres for all time steps of
    the current simulation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass to upscale magnitudes and with all passed
        # parameters.
        super().__init__(
            magnitude_factor=_MAGNITUDE_FACTOR_CURRENT,
            *args, **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property_cached
    def x(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the X
          component of the intracellular current density vector spatially
          situated at the center of each grid space for this time step.
        '''

        return self._get_component(self._sim.I_cell_x_time)


    @property_cached
    def y(self) -> ndarray:
        '''
        Two-dimensional Numpy array whose:

        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the Y
          component of the intracellular current density vector spatially
          situated at the center of each grid space for this time step.
        '''

        return self._get_component(self._sim.I_cell_y_time)

    # ..................{ GETTERS                            }..................
    @type_check
    def _get_component(self, times_currents_centered: SequenceTypes) -> ndarray:
        '''
        Interpolate the passed array of the X and Y components of all
        intracellular current density vectors in this vector field for all time
        steps from heterogenous cell centers onto the homogenous world grid.

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
            * First dimension indexes each simulation time step.
            * Second dimension indexes the X or Y component of each
              intracellular current density vector computed at the center of
              each cell for this time step.

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
