#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all timed event classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractmethod
from betse.util.type.types import type_check, NumericSimpleTypes

# ....................{ BASE                               }....................
class SimEventABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all timed event classes.

    Instances of this class parameterize an exogenous event (e.g., cell removal,
    voltage change) to be triggered at some time step(s) of the simulation.

    Attributes
    ----------
    _is_fired : bool
        ``True`` only if this event's :meth:`fire` method has already been
        called for a previous time step of the appropriate simulation phase.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this event.
        '''

        # Classify all passed parameters.
        self._is_fired = False

    # ..................{ PROPERTIES                         }..................
    # Read-only properties, preventing callers from setting these attributes.

    @property
    def is_fired(self) -> bool:
        '''
        ``True`` only if this event's :meth:`fire` method has already been
        called for a previous time step of the appropriate simulation phase.
        '''

        return self._is_fired

    # ..................{ FIRERS                             }..................
    @type_check
    def fire(self) -> None:
        '''
        Apply this event.

        This method notes this event to have been applied, ensuring the
        :meth:`is_fired` property to subsequently report ``True``.
        '''

        self._is_fired = True

# ....................{ SUBCLASSES                         }....................
class SimEventSpikeABC(SimEventABC):
    '''
    Abstract base class of all classes describing simulation events occurring at
    only a single time step (rather than over a range of time steps).

    Attributes
    ----------
    _time_step : NumericSimpleTypes
        Time step in seconds (s) at which to trigger this action.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, time_step: NumericSimpleTypes) -> None:
        '''
        Initialize this event.

        Parameters
        ----------
        time_step : NumericSimpleTypes
            Time step in seconds (s) at which to trigger this action.
        '''

        # Initialize our superclass.
        super().__init__()

        # Classify all passed parameters.
        self._time_step = time_step


class SimEventPulseABC(SimEventABC):
    '''
    Abstract base class of all classes describing simulation events occurring
    over a range of time steps (rather than at only a single time step).

    Attributes
    ----------
    start_time : NumericSimpleTypes
        Time step (s) at which to begin triggering this event.
    stop_time : NumericSimpleTypes
        Time step (s) at which to cease triggering this event.
    step_rate : NumericSimpleTypes
        Slope of the pair of step functions guaranteeing smooth continuity
        between the background function and this event. Each step function is
        the mirror image of the other reflected across the Y axis. These are:
        * A "step" up from the background function to this event, whose
          mid-point is centered at ``start_time``.
        * A "step" down from this event back to the background function, whose
          mid-point is centered at ``stop_time``.
        If the background function is time-dependent, this slope is a **rate**
        (i.e., change over time). For the :class:`SimEventPulseVoltage` subclass, for
        example, this is the rate in voltage per seconds (V/s) at which:
        * The background voltage is first increased to the peak voltage.
        * The peak voltage is later decreased to the background voltage.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        start_time: NumericSimpleTypes,
        stop_time: NumericSimpleTypes,
        step_rate: NumericSimpleTypes,
    ) -> None:

        # Initialize our superclass.
        super().__init__()

        # Classify all passed parameters.
        self.start_time = start_time
        self.stop_time = stop_time
        self.step_rate = step_rate
