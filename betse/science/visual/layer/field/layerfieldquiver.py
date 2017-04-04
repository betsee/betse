#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying vector components as quiver plots onto the
cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.export import expmath
from betse.science.visual.layer.field.layerfieldabc import (
    LayerCellsFieldColorlessABC)
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class LayerCellsFieldQuiver(LayerCellsFieldColorlessABC):
    '''
    Layer subclass plotting the most significant X and Y components of a single
    vector field (e.g., electric field) onto the cell cluster for one on more
    simulation time steps.

    Attributes
    ----------
    _quiver_plot : matplotlib.quiver.Quiver
        Ouiver plot of all vector components previously plotted for the prior
        time step if any or `None` otherwise.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass.
        super().__init__(*args, **kwargs)

        # Default all remaining instance variables.
        self._quiver_plot = None

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:
        '''
        Layer the most significant X and Y components of this vector field for
        the first time step onto the figure axes of the current plot or
        animation.
        '''

        # Vector field whose X and Y components are spatially situated at cell
        # centres.
        field = self._field.times_cells_centre

        # Arrays of the upscaled X and Y coordinates of all cell centres.
        cells_centre_x = expmath.upscale_cell_coordinates(
            self._phase.cells.cell_centres[:,0])
        cells_centre_y = expmath.upscale_cell_coordinates(
            self._phase.cells.cell_centres[:,1])

        # Ouiver plot of all vector components plotted for this time step. See
        # the matplotlib.quiver.quiver() docstring for further details.
        self._quiver_plot = self._visual.axes.quiver(
            # Positional arguments. Thanks to internal flaws in the
            # matplotlib.quiver._parse_args() function parsing arguments passed
            # to the matplotlib.axes.quiver() method called here, the first four
            # arguments *MUST* be passed as positional arguments.

            # X and Y coordinates of all cell centres.
            cells_centre_x,
            cells_centre_y,

            # Normalized X and Y vector field components for this time step.
            field.unit_x[self._visual.time_step],
            field.unit_y[self._visual.time_step],

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

            #FIXME: This appears to be ignored. Is this still required?
            # zorder=10,
        )


    def _layer_next(self) -> None:
        '''
        Layer the most significant X and Y components of this vector field for
        the next time step onto the figure axes of the current plot or
        animation.
        '''

        # Vector field whose X and Y components are spatially situated at cell
        # centres.
        field = self._field.times_cells_centre

        # Arrays of all normalized X and Y components of this vector field for
        # this time step.
        field_unit_x = field.unit_x[self._visual.time_step]
        field_unit_y = field.unit_y[self._visual.time_step]

        # Replace all normalized X and Y components previously plotted for the
        # prior time step by these components.
        self._quiver_plot.set_UVC(U=field_unit_x, V=field_unit_y)
