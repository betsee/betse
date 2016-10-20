#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses spatially overlaying streamlines onto the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import abstractmethod
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.types import type_check

# ....................{ CLASSES ~ base                     }....................
class LayerCellsStreamCurrent(LayerCellsABC):
    '''
    Abstract base class of all layer subclasses plotting streamlines of
    electrical current density onto the cell cluster.

    Such layers are somewhat more computationally expensive in both space and
    time than the average layer. For each plot or animation frame to be layered
    with streamlines, the subclass solves an internal fluid simulation of the
    current density through this cell cluster specific to this frame.

    Attributes
    ----------
    _time_current_x : SequenceTypes
        Two-dimensional sequence whose:
        * First dimension indexes all time steps of this simulation.
        * Second dimension indexes the X components of all current density
          vectors computed for the corresponding time step.
    _time_current_y : SequenceTypes
        Two-dimensional sequence of the same format as `time_current_x`,
        replacing "X components" by "Y components".
    _time_current_magnitude : SequenceTypes
        Two-dimensional sequence whose:
        * First dimension indexes all time steps of this simulation.
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

        # For efficiency, coerce these sequences into arrays.
        self._time_currents_x = np.asarray(self._time_current_x)
        self._time_currents_y = np.asarray(self._time_current_y)

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
    @type_check
    def layer_first(self) -> None:
        '''
        Layer the current density for the current current time step and cell
        cluster onto the figure axes of the passed plot or animation as an
        arbitrary number of streamlines.
        '''

        # Current density magnitudes for this time step.
        currents_magnitude = self._time_currents_magnitude[
            self._visual.time_step]

        #FIXME: Implement us up, please.

        # Erase the prior frame's overlay and streamplot this frame's overlay.
        # self._current_density_stream_plot = self._plot_stream(
        #     old_stream_plot=self._current_density_stream_plot,
        #     x=self._current_density_x_time_series[frame_number] / Jmag_M,
        #     y=self._current_density_y_time_series[frame_number] / Jmag_M,
        #     magnitude=Jmag_M,
        # )

# ....................{ CLASSES ~ sub                      }....................
class LayerCellsStreamCurrentIntraExtra(LayerCellsStreamCurrent):
    '''
    Layer subclass plotting streamlines of the current density of all
    intracellular and extracellular spaces onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _prep_time_currents_components(self) -> None:

        self._time_currents_x = self._visual.sim.I_tot_x_time
        self._time_currents_y = self._visual.sim.I_tot_y_time


class LayerCellsStreamCurrentIntra(LayerCellsStreamCurrent):
    '''
    Layer subclass plotting streamlines of the current density of all
    intracellular spaces (e.g., gap junctions) onto the cell cluster.
    '''

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Implement us up, please.
    def _prep_time_currents_components(self) -> None:

        pass
