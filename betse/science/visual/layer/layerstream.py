#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying streamlines onto the current cell cluster.
'''

#FIXME: Optimize. The current _get_time_currents() approach is extremely
#inefficient, as this entire vector field will be recomputed for each plot and
#animation when current overlays are enabled. The simpliest alternative is to
#add two new properties to the "Simulator" class:
#
#* Simulator.time_currents_intra(), returning the vector field...
#FIXME: Oh, wait. These vector fields require the "sim", "cells", and "p"
#parameters -- which, due to arguably poor Simulator design, are currently
#unavailable as attributes there. Very well. We see little choice but to finally
#define a pipeline class. As these vector fields are currently only required by
#plots and animations, sequestering these fields to that pipeline class would be
#trivial. However, we would then need to pass the same instance of that class to
#*EVERY* instantiated plot and animation.
#
#Doing so is ultimately trivial but tedious and hence deferred to another day.

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseVectorException
from betse.science.vector.field.fieldabc import VectorFieldABC
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.types import type_check
from matplotlib.patches import FancyArrowPatch

# ....................{ SUBCLASSES                         }....................
#FIXME: Generalize this layer to accept a vector field whose X and Y components
#are *NOT* spatially situated at square grid spaces, presumably by interpolating
#from the coordinate system of these components onto square grid spaces.
#FIXME: To do so, subclass from the new "LayerCellsVectorFieldABC" superclass
#instead.

class LayerCellsStream(LayerCellsABC):
    '''
    Layer subclass plotting streamlines of a single vector field whose X and Y
    components are spatially situated at square grid spaces (e.g.,
    intracellular current density) for one on more simulation time steps.

    This layer is somewhat more computationally expensive in both space and time
    than the average layer. For each plot or animation frame to be layered with
    streamlines, an internal fluid simulation of the density of the desired
    vector field through the cell cluster specific to this frame is solved.

    Attributes
    ----------
    _field : VectorFieldABC
        Vector field of all velocity vectors spatially situated at square grid
        spaces for all simulation time steps to be streamplotted by this layer.
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


    def prep(self, *args, **kwargs) -> None:
        '''
        Prepare this layer to be layered onto the passed plot or animation.

        Parameters
        ----------
        All parameters are passed to the :meth:`LayerCellsABC.__init__` method
        as is.

        Raises
        ----------
        BetseVectorException
            If the X and Y components of this vector field are *not* spatially
            situated at square grid spaces.
        '''

        # Prepare our superclass with the passed parameters.
        super().prep(*args, **kwargs)

        # If the X and Y components of this vector field are *not* spatially
        # situated at square grid spaces, raise an exception. Ideally, this
        # validation would be performed in the constructor. Since the
        # "self._visual" attribute is only defined by the prior call to our
        # superclass implementation, this validation is delayed to here.
        if not (
            len(self._field.x[0]) == len(self._visual.cells.X) and
            len(self._field.y[0]) == len(self._visual.cells.Y)
        ):
            raise BetseVectorException(
                'Vector field not defined at square grid spaces (i.e., '
                'expected dimensions {} x {} but '
                'received dimensions {} x {}).'.format(
                    len(self._visual.cells.X),
                    len(self._visual.cells.Y),
                    len(self._field.x[0]), len(self._field.y[0])))

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
