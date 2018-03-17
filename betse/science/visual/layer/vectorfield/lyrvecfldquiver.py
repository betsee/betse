#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying vector components as quiver plots onto the
cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.visual.layer.vectorfield.lyrvecfldabc import (
    LayerCellsFieldColorlessABC)
from betse.util.type.decorator.deccls import abstractproperty
# from betse.util.type.types import type_check
from numpy import ndarray

# ....................{ SUPERCLASSES                       }....................
class LayerCellsFieldQuiverABC(LayerCellsFieldColorlessABC):
    '''
    Layer abstract subclass plotting the most significant X and Y components of
    a single vector field spatially situated according to this layer subclass
    onto the cell cluster for one on more time steps.

    Attributes
    ----------
    _quiver_plot : matplotlib.quiver.Quiver
        Ouiver plot of all vector components previously plotted for the prior
        time step if any *or* ``None`` otherwise.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass.
        super().__init__(*args, **kwargs)

        # Default all remaining instance variables.
        self._quiver_plot = None

    # ..................{ SUPERCLASS                         }..................
    def _layer_next(self) -> None:

        # Replace all X and Y components of this vector field plotted for the
        # prior time step by those for this time step.
        self._quiver_plot.set_UVC(U=self._field_time_x, V=self._field_time_y)

    # ..................{ SUBCLASS                           }..................
    @abstractproperty
    def _field_time_x(self) -> ndarray:
        '''
        One-dimensional Numpy array of all X components of this vector field
        spatially situated spatially situated according to this layer subclass
        for the current time step.
        '''

        pass


    @abstractproperty
    def _field_time_y(self) -> ndarray:
        '''
        One-dimensional Numpy array of all Y components of this vector field
        spatially situated spatially situated according to this layer subclass
        for the current time step.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class LayerCellsFieldQuiverCells(LayerCellsFieldQuiverABC):
    '''
    Layer subclass plotting the most significant X and Y components of a single
    vector field spatially situated at cell centres (e.g., intracellular
    electric field) onto the cell cluster for one on more time steps.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:

        # Ouiver plot of all vector components plotted for this time step. See
        # the matplotlib.quiver.quiver() docstring for further details.
        self._quiver_plot = self._visual.axes.quiver(
            # Positional arguments. Thanks to internal flaws in the
            # matplotlib.quiver._parse_args() function parsing arguments passed
            # to the matplotlib.axes.quiver() method called here, the first four
            # arguments *MUST* be passed as positional arguments.

            # Upscaled X and Y coordinates of all cell centres.
            self._phase.cache.upscaled.cells_centre_x,
            self._phase.cache.upscaled.cells_centre_y,

            # Normalized X and Y components of this vector field spatially
            # situated at cell centres for this time step.
            self._field_time_x,
            self._field_time_y,

            # Keyword arguments. All remaining arguments *MUST* be passed as
            # keyword arguments.

            # Matplotlib-specific color code of all vector arrows.
            color=self._phase.p.vcolor,

            # Multiples of the width and height (respectively) of vector arrow
            # shafts by which to scale the width and height of vector arrow
            # heads. These settings default to 3 and 5 (respectively).
            headwidth=5,
            headlength=7,

            # The portion of each vector arrow to situate at the X and Y
            # coordinates of the corresponding cell centre.
            pivot='middle',

            # Scale vector arrows such that arrow size increases as the
            # user-defined zoom level increases in either X or Y dimensions.
            # units='xy',
            units='x',

            # Z-order of this plot with respect to other artists.
            zorder=self._zorder,
        )

    # ..................{ SUPERCLASS ~ property              }..................
    @property
    def _field_time_x(self) -> ndarray:
        return self._field.times_cells_centre.unit_x[self._visual.time_step]

    @property
    def _field_time_y(self) -> ndarray:
        return self._field.times_cells_centre.unit_y[self._visual.time_step]


class LayerCellsFieldQuiverGrids(LayerCellsFieldQuiverABC):
    '''
    Layer subclass plotting the most significant X and Y components of a single
    vector field spatially situated at environmental grid space centres (e.g.,
    extracellular electric field) onto the cell cluster for one on more time
    steps.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:

        # Ouiver plot of all vector components plotted for this time step. See
        # the matplotlib.quiver.quiver() docstring for further details.
        self._quiver_plot = self._visual.axes.quiver(
            # Positional arguments. See
            # LayerCellsFieldQuiverCells._layer_first() for further discussion.

            # Upscaled X and Y coordinates of all grid space centres.
            self._phase.cache.upscaled.grids_centre_x,
            self._phase.cache.upscaled.grids_centre_y,

            # Normalized X and Y components of this vector field spatially
            # situated at grid space centres for this time step.
            self._field_time_x,
            self._field_time_y,

            #FIXME: DRY. All following keyword arguments are duplicated from the
            #LayerCellsFieldQuiverCells._layer_first() method, which is bad. In
            #fact, the only difference between this and that class is the use
            #of "grids_centre"-centric arrays above and below. So, four
            #differences in total. Can we generalize this further?

            # Keyword arguments. All remaining arguments *MUST* be passed as
            # keyword arguments.

            # Matplotlib-specific color code of all vector arrows.
            color=self._phase.p.vcolor,

            # Multiples of the width and height (respectively) of vector arrow
            # shafts by which to scale the width and height of vector arrow
            # heads. These settings default to 3 and 5 (respectively).
            headwidth=5,
            headlength=7,

            # The portion of each vector arrow to situate at the X and Y
            # coordinates of the corresponding cell centre.
            pivot='middle',

            # Scale vector arrows such that arrow size increases as the
            # user-defined zoom level increases in either X or Y dimensions.
            # units='xy',
            units='x',

            # Z-order of this plot with respect to other artists.
            zorder=self._zorder,
        )

    # ..................{ SUPERCLASS ~ property              }..................
    @property
    def _field_time_x(self) -> ndarray:

        field = self._field.times_grids_centre
        return (
            field.x[self._visual.time_step] /
            field.magnitudes[self._visual.time_step].max())


    @property
    def _field_time_y(self) -> ndarray:

        field = self._field.times_grids_centre
        return (
            field.y[self._visual.time_step] /
            field.magnitudes[self._visual.time_step].max())


class LayerCellsFieldQuiverMembranes(LayerCellsFieldQuiverABC):
    '''
    Layer subclass plotting the most significant X and Y components of a single
    vector field spatially situated at cell membrane midpoints (e.g.,
    microtubules) onto the cell cluster for one on more time steps.

    Attributes
    ----------
    _cells_radius : ndarray
        One-dimensional Numpy array indexing each cell such that each element is
        the upscaled radius of that cell.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass.
        super().__init__(*args, **kwargs)

        # Default all remaining instance variables.
        self._cells_radius = None


    def prep(self, *args, **kwargs) -> None:

        # Prepare our superclass.
        super().prep(*args, **kwargs)

        # Upscale all cell radii.
        self._cells_radius = expmath.upscale_coordinates(
            self._phase.cells.R[self._phase.cells.mem_to_cells])

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:

        # Localize attributes for brevity.
        cells = self._phase.cells

        # Ouiver plot of all vector components plotted for this time step. See
        # the matplotlib.quiver.quiver() docstring for further details.
        self._quiver_plot = self._visual.axes.quiver(
            # Positional arguments. See
            # LayerCellsFieldQuiverCells._layer_first() for further discussion.

            # Upscaled X and Y coordinates of all cell membrane midpoints.
            expmath.upscale_coordinates(
                cells.cell_centres[:, 0][cells.mem_to_cells]),
            expmath.upscale_coordinates(
                cells.cell_centres[:, 1][cells.mem_to_cells]),

            # X and Y components of this vector field spatially situated at cell
            # membrane midpoints extending to cell centres for this time step.
            self._field_time_x,
            self._field_time_y,

            # Keyword arguments. All remaining arguments *MUST* be passed as
            # keyword arguments.

            # Matplotlib-specific color code of all vector arrows.
            color=self._phase.p.vcolor,

            # Number of data units per arrow length unit.
            scale=expmath.upscale_coordinates(self._phase.p.wsx * 0.8),

            # Z-order of this plot with respect to other artists.
            zorder=self._zorder,
        )

    # ..................{ SUPERCLASS ~ property              }..................
    @property
    def _field_time_x(self) -> ndarray:
        '''
        X components of this vector field spatially situated at cell membrane
        midpoints extending exactly to cell centres for the current time step,
        calculated by first normalizing and then extending each such component
        by the upscaled radius of the cell containing this membrane.
        '''

        return self._field.times_membranes_midpoint.unit_x[
            self._visual.time_step] * self._cells_radius


    @property
    def _field_time_y(self) -> ndarray:
        '''
        Y components of this vector field spatially situated at cell membrane
        midpoints extending exactly to cell centres for the current time step.
        '''

        return self._field.times_membranes_midpoint.unit_y[
            self._visual.time_step] * self._cells_radius
