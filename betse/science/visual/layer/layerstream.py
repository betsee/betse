#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying streamlines onto the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import abstractmethod
from betse.exceptions import BetseSimConfigException
from betse.lib.numpy import arrays
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type import sequences
from betse.util.type.types import type_check, SequenceTypes
from matplotlib.patches import FancyArrowPatch
from scipy import interpolate

# ....................{ SUPERCLASSES                       }....................
class LayerCellsStream(LayerCellsABC):
    '''
    Abstract base class of all layer subclasses plotting streamlines of a single
    modelled vector field (e.g., intracellular current) for one on more
    simulation time steps.

    Such layers are somewhat more computationally expensive in both space and
    time than the average layer. For each plot or animation frame to be layered
    with streamlines, an internal fluid simulation of the density of the desired
    vector field through the cell cluster specific to this frame is solved.

    Attributes
    ----------
    _streamplot : StreamplotSet
        Streamplot of all streamlines previously plotted for the prior time step
        if any or `None` otherwise, temporarily preserved for only one time step
        to permit its removal prior to plotting a new streamplot for the current
        time step.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._streamplot = None

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def _get_velocities_x(self) -> SequenceTypes:
        '''
        Array of the X components of all velocity vectors in this vector field
        for the current time step.
        '''

        pass


    @abstractmethod
    def _get_velocities_y(self) -> SequenceTypes:
        '''
        Array of the Y components of all velocity vectors in this vector field
        for the current time step.
        '''

        pass


    @abstractmethod
    def _get_velocities_magnitude(self) -> SequenceTypes:
        '''
        Array of the magnitudes of all velocity vectors in this vector field for
        the current time step.
        '''

        pass

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:
        '''
        Simulate and layer streamlines of a single modelled vector field (e.g.,
        intracellular current) for the next time step onto the figure axes of
        the current plot or animation.
        '''

        # Arrays of the upscaled X and Y coordinates of all grid points.
        grid_x = visuals.upscale_cell_coordinates(self._visual.cells.X)
        grid_y = visuals.upscale_cell_coordinates(self._visual.cells.Y)

        # Arrays of all X and Y components and magnitudes of this vector field
        # for this time step.
        velocities_x = self._get_velocities_x()
        velocities_y = self._get_velocities_y()
        velocities_magnitude = self._get_velocities_magnitude()

        # Arrays of all normalized X and Y components of this vector field for
        # this time step.
        normalized_velocities_x = velocities_x / velocities_magnitude
        normalized_velocities_y = velocities_y / velocities_magnitude

        # Maximum magnitude of this vector field for this time step.
        velocities_magnitude_max = np.max(velocities_magnitude)

        # Array of display-specific line widths for all vectors of this vector
        # field for this time step.
        streamlines_width = (
            3.0 * velocities_magnitude / velocities_magnitude_max) + 0.5

        # Streamplot of all streamlines plotted for this time step. See the
        # matplotlib.streamplot.streamplot() docstring for further details.
        self._streamplot = self._visual.axes.streamplot(
            # X and Y coordinates of all grid points.
            x=grid_x,
            y=grid_y,

            # X and Y coomponents of all vector field velocities.
            u=normalized_velocities_x,
            v=normalized_velocities_y,

            # Matplotlib-specific color code of all streamlines.
            color=self._visual.p.vcolor,

            # Density of streamlines in both the X and Y dimensions.
            density=self._visual.p.stream_density,

            # Line widths of all streamlines.
            linewidth=streamlines_width,

            # Factor by which to upscale the size of all streamline arrowheads.
            arrowsize=1.5,
        )


    def _layer_next(self) -> None:
        '''
        Simulate and layer streamlines of a single modelled vector field (e.g.,
        intracellular current) for the first time step onto the figure axes of
        the current plot or animation.
        '''

        # Remove all streamlines plotted for the prior time step.
        self._streamplot.lines.remove()

        # If this Matplotlib version supports removing the set of all streamline
        # arrowheads plotted for the prior time step, do so.
        try:
            self._streamplot.arrows.remove()
        # Else, these arrowheads *MUST* be manually erased by iterating over all
        # patch objects and preserving all non-arrowhead patches. Doing so also
        # removes all arrowhead patches of other streamplots already plotted for
        # this time step and is hence non-ideal. But no alternatives exist.
        except NotImplementedError:
            self._visual.axes.patches = [
                patch
                for patch in self._visual.axes.patches
                if not isinstance(patch, FancyArrowPatch)
            ]

        # Replot this streamplot for this time step.
        self._layer_first()


class LayerCellsStreamCurrent(LayerCellsStream):
    '''
    Abstract base class of all layer subclasses plotting streamlines of
    electrical current density onto the cell cluster.

    Such layers are somewhat more computationally expensive in both space and
    time than the average layer. For each plot or animation frame to be layered
    with streamlines, the subclass solves an internal fluid simulation of the
    current density through this cell cluster specific to this frame.

    Attributes
    ----------
    _time_currents_x : SequenceTypes
        Two-dimensional sequence whose:
        * First dimension indexes each time step of this simulation.
        * Second dimension indexes the X components of all current density
          vectors computed for the corresponding time step.
    _time_currents_y : SequenceTypes
        Two-dimensional sequence of the same format as `time_current_x`,
        replacing "X components" by "Y components".
    _time_currents_magnitude : SequenceTypes
        Two-dimensional sequence whose:
        * First dimension indexes each time step of this simulation.
        * Second dimension indexes the magnitudes (commonly referred to as
          `Jmag_M` elsewhere in the codebase) of all current density vectors
          computed for the corresponding time step, in units of uA/cm2.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._time_currents_x = None
        self._time_currents_y = None
        self._time_currents_magnitude = None


    def prep(self, *args, **kwargs) -> None:

        # Prepare our superclass with all passed parameters.
        super().prep(*args, **kwargs)

        # Define the "_time_currents_x" and "_time_currents_y" attributes.
        self._prep_time_currents_components()

        # For efficiency, coerce these possibly non-array sequences into arrays.
        self._time_currents_x = arrays.from_sequence(self._time_currents_x)
        self._time_currents_y = arrays.from_sequence(self._time_currents_y)

        # If either of these arrays is empty (e.g., due to erroneously
        # attempting to layer extracellular current with extracellular spaces
        # disabled), raise an exception.
        sequences.die_if_empty(
            self._time_currents_x, label='X current components')
        sequences.die_if_empty(
            self._time_currents_y, label='Y current components')

        # Array of current density magnitude computed from these arrays:
        #
        # * Incremented by a negligible positive value approximately equal to 0,
        #   avoiding division-by-zero errors when the layer() method
        #   subsequently divides by elements of this array.
        # * Multiplied by 100, yielding magnitude in units of uA/cm2.
        self._time_currents_magnitude = 1e-15 + 100 * np.sqrt(
            self._time_currents_x ** 2 +
            self._time_currents_y ** 2)

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def _prep_time_currents_components(self) -> None:
        '''
        Define the :attr:`_time_currents_x` and :attr:`_time_currents_y` arrays.
        '''

        pass

    # ..................{ SUPERCLASS                         }..................
    def _get_velocities_x(self) -> SequenceTypes:
        '''
        Array of the X components of all velocity vectors in this vector field
        for the current time step.
        '''

        return self._time_currents_x[self._visual.time_step]


    def _get_velocities_y(self) -> SequenceTypes:
        '''
        Array of the Y components of all velocity vectors in this vector field
        for the current time step.
        '''

        return self._time_currents_y[self._visual.time_step]


    def _get_velocities_magnitude(self) -> SequenceTypes:
        '''
        Array of the magnitudes of all velocity vectors in this vector field for
        the current time step.
        '''

        return self._time_currents_magnitude[self._visual.time_step]

# ....................{ SUBCLASSES                         }....................
class LayerCellsStreamCurrentIntraExtra(LayerCellsStreamCurrent):
    '''
    Layer subclass plotting streamlines of the current density of all
    intracellular and extracellular spaces onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _prep_time_currents_components(self) -> None:
        '''
        Define the :attr:`_time_currents_x` and :attr:`_time_currents_y` arrays
        to simply be synonyms of the :attr:`Simulator.sim.I_tot_x_time` and
        :attr:`Simulator.sim.I_tot_y_time` arrays (_respectively_).

        Raises
        ----------
        BetseSimConfigException
            If extracellular spaces are disabled.
        '''

        # If extracellular spaces are disabled, raise an exception.
        if not self._visual.p.sim_ECM:
            raise BetseSimConfigException('Extracellular spaces disabled.')

        self._time_currents_x = self._visual.sim.I_tot_x_time
        self._time_currents_y = self._visual.sim.I_tot_y_time


class LayerCellsStreamCurrentIntra(LayerCellsStreamCurrent):
    '''
    Layer subclass plotting streamlines of the current density of all
    intracellular spaces (e.g., gap junctions) onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _prep_time_currents_components(self) -> None:
        '''
        Interpolate the :attr:`_time_currents_x` and :attr:`_time_currents_y`
        arrays from the current density of only intracellular spaces defined
        at cell centers to the X/Y grid of this cell cluster's environment.

        In the case of the :class:`LayerCellsStreamCurrentIntraExtra` subclass,
        the current density of all intracellular and extracellular spaces is
        precomputed by the :class:`Simulator` class and hence reusable as is.
        In the case of this subclass, however, the current density of only
        intracellular spaces is _not_ precomputed elsewhere and hence must be
        computed here.
        '''

        # print('!!!!in Intra')
        # Two-dimensional tuple of the X and Y components of all cell centers.
        dimensions_cells_center = (
            self._visual.cells.cell_centres[:, 0],
            self._visual.cells.cell_centres[:, 1])

        # Two-dimensional tuple of the X and Y components of all grid points.
        dimensions_grid_points = (self._visual.cells.X, self._visual.cells.Y)

        # Two-dimensional lists of the X and Y components of all current density
        # vectors over all time steps.
        self._time_currents_x = self._prep_currents_component(
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
    def _prep_currents_component(
        self,
        times_currents_center: SequenceTypes,
        dimensions_cells_center: SequenceTypes,
        dimensions_grid_points: SequenceTypes,
    ) -> SequenceTypes:
        '''
        Interpolate the passed array of the X or Y components of all
        intracellular current density vectors defined at cell centers onto the
        passed X/Y grid of this cell cluster's environment.

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
            currents_grid = self._visual.cells.maskECM * interpolate.griddata(
                # Tuple of X and Y coordinates to interpolate from.
                points=dimensions_cells_center,

                # Tuple of X and Y coordinates to interpolate onto.
                xi=dimensions_grid_points,

                # Array of X or Y current vector components to interpolate.
                values=currents_center,

                # Machine-readable string specifying the type of interpolation
                # to perform, established by the current configuration.
                method=self._visual.p.interp_type,

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
