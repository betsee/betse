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
from abc import abstractmethod, abstractproperty

import numpy as np
from betse.science.vector.fieldelectric import (
    VectorFieldSimmedABC,
    VectorFieldCurrentIntra,
    VectorFieldCurrentIntraExtra,
)
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.callables import property_cached
from betse.util.type.types import type_check, SequenceTypes
from matplotlib.patches import FancyArrowPatch


# ....................{ SUPERCLASSES                       }....................
class LayerCellsStreamABC(LayerCellsABC):
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

        # Initialize our superclass.
        super().__init__()

        # Default instance attributes.
        self._streamplot = None

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def _get_velocities_x(self) -> SequenceTypes:
        '''
        Numpy array of the X components of all velocity vectors in this vector
        field for the current time step.
        '''

        pass


    @abstractmethod
    def _get_velocities_y(self) -> SequenceTypes:
        '''
        Numpy array of the Y components of all velocity vectors in this vector
        field for the current time step.
        '''

        pass


    @abstractmethod
    def _get_velocities_magnitudes(self) -> SequenceTypes:
        '''
        Numpy array of the magnitudes of all velocity vectors in this vector
        field for the current time step.
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
        velocities_magnitudes = self._get_velocities_magnitudes()

        # Arrays of all normalized X and Y components of this vector field for
        # this time step.
        normalized_velocities_x = velocities_x / velocities_magnitudes
        normalized_velocities_y = velocities_y / velocities_magnitudes

        # Maximum magnitude of this vector field for this time step.
        velocities_magnitude_max = np.max(velocities_magnitudes)

        # Array of display-specific line widths for all vectors of this vector
        # field for this time step.
        streamlines_width = (
            3.0 * velocities_magnitudes / velocities_magnitude_max) + 0.5

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

            # FIXME: for stillframes an arrow size of about 5.0 is best; for rendered video these
            # blow up and a size of 3.0 is better. Not that this is important enough to ever get to it :)
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


class LayerCellsStreamCurrentABC(LayerCellsStreamABC):
    '''
    Abstract base class of all layer subclasses plotting streamlines of
    electrical current density onto the cell cluster.

    Such layers are somewhat more computationally expensive in both space and
    time than the average layer. For each plot or animation frame to be layered
    with streamlines, the subclass solves an internal fluid simulation of the
    current density through this cell cluster specific to this frame.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _get_velocities_x(self) -> SequenceTypes:
        '''
        Numpy array of the X components of all velocity vectors in this vector
        field for the current time step.
        '''

        return self._time_currents.x[self._visual.time_step]


    def _get_velocities_y(self) -> SequenceTypes:
        '''
        Numpy array of the Y components of all velocity vectors in this vector
        field for the current time step.
        '''

        return self._time_currents.y[self._visual.time_step]


    def _get_velocities_magnitudes(self) -> SequenceTypes:
        '''
        Numpy array of the magnitudes of all velocity vectors in this vector
        field for the current time step.
        '''

        return self._time_currents.magnitudes[self._visual.time_step]

    # ..................{ SUBCLASS                           }..................
    @abstractproperty
    def _time_currents(self) -> VectorFieldSimmedABC:
        '''
        Vector field of the current densities of all intracellular and/or
        extracellular spaces spatially situated at grid space centres for all
        time steps of the current simulation.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class LayerCellsStreamCurrentIntraExtra(LayerCellsStreamCurrentABC):
    '''
    Layer plotting streamlines of the current density of all intracellular and
    extracellular spaces onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    @property_cached
    def _time_currents(self) -> VectorFieldSimmedABC:
        return VectorFieldCurrentIntraExtra(
            sim=self._visual._sim,
            cells=self._visual._cells,
            p=self._visual._p,
        )


class LayerCellsStreamCurrentIntra(LayerCellsStreamCurrentABC):
    '''
    Layer plotting streamlines of the current density of only all intracellular
    spaces (e.g., gap junctions) onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    @property_cached
    def _time_currents(self) -> VectorFieldSimmedABC:
        return VectorFieldCurrentIntra(
            sim=self._visual._sim,
            cells=self._visual._cells,
            p=self._visual._p,
        )
