#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying vector components onto the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.vector.field.fieldabc import VectorFieldABC
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.types import type_check
from matplotlib.patches import FancyArrowPatch

# ....................{ SUBCLASSES                         }....................
class LayerCellsQuiver(LayerCellsABC):
    '''
    Layer subclass plotting vector components of a single vector field (e.g.,
    intracellular electric field) for one on more simulation time steps.

    This layer is somewhat more computationally expensive in both space and time
    than the average layer. For each plot or animation frame to be layered with
    streamlines, an internal fluid simulation of the density of the desired
    vector field through the cell cluster specific to this frame is solved.

    Attributes
    ----------
    _field : VectorFieldABC
        Vector field of all velocity vectors for all simulation time steps to be
        streamplotted by this layer.
    _streamplot : StreamplotSet
        Streamplot of all streamlines previously plotted for the prior time step
        if any or `None` otherwise, temporarily preserved for only one time step
        to permit its removal prior to plotting a new streamplot for the current
        time step.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, field: VectorFieldABC) -> None:
        '''
        Initialize this layer.

        Parameters
        ----------
        field : VectorFieldABC
            Vector field of all velocity vectors for all simulation time steps
            to be streamplotted by this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify the passed parameter.
        self._field = field

        # Default all remaining instance variables.
        self._streamplot = None

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

        # Arrays of all magnitudes *AND* normalized X and Y components of this
        # vector field for this time step.
        field_magnitudes = self._field.magnitudes[self._visual.time_step]
        field_unit_x = self._field.unit_x[self._visual.time_step]
        field_unit_y = self._field.unit_y[self._visual.time_step]

        # Maximum magnitude of this vector field for this time step.
        field_magnitude_max = np.max(field_magnitudes)

        # One-dimensional array of the visual widths of all vectors of this
        # vector field for this time step.
        streamlines_width = (
            3.0 * field_magnitudes / field_magnitude_max) + 0.5

        # Streamplot of all streamlines plotted for this time step. See the
        # matplotlib.streamplot.streamplot() docstring for further details.
        self._streamplot = self._visual.axes.streamplot(
            # X and Y coordinates of all grid points.
            x=grid_x,
            y=grid_y,

            # X and Y normalized coomponents of all vector field velocities.
            u=field_unit_x,
            v=field_unit_y,

            # Matplotlib-specific color code of all streamlines.
            color=self._visual.p.vcolor,

            # Density of streamlines in both the X and Y dimensions.
            density=self._visual.p.stream_density,

            # Line widths of all streamlines.
            linewidth=streamlines_width,

            #FIXME: For still frames, an arrow size of about 5.0 is best; for
            #rendered video, these blow up and a size of 3.0 is better. Not that
            #this is important enough to ever get to it. :)

            # Factor by which to upscale the size of all streamline arrowheads.
            arrowsize=3.0,
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
