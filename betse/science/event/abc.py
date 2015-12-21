#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes for timed event classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractstaticmethod
from betse.util.type import types

# ....................{ BASE                               }....................
class Event(object, metaclass = ABCMeta):
    '''
    Abstract base class of all timed event classes.

    Instances of this class parameterize an exogenous event (e.g., cell removal,
    voltage change) to be triggered at some time step(s) of the simulation.

    Attributes
    ----------------------------
    name : str
        Unique name of the current profile.
    '''

    # ..................{ ABSTRACT ~ static                  }..................
    @abstractstaticmethod
    def make(params: 'Parameters') -> 'Event':
        '''
        Factory method producing a concrete instance of this abstract base class
        from the passed tissue simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
             Current tissue simulation configuration.

        Returns
        ----------------------------
        Event
            Concrete instance of this abstract base class.
        '''
        pass

# ....................{ PERIOD                             }....................
class EventPeriod(Event):
    '''
    Abstract base class of all classes describing simulation events occurring
    over a period (rather than single point) of time.

    Attributes
    ----------------------------
    start_time : float
        Time point in seconds at which to begin triggering this event.
    stop_time : float
        Time point in seconds at which to cease triggering this event.
    step_width : float
        Time period in seconds of the width of the step function providing
        smooth continuity between the underlying time-dependent function and
        the event overlayed onto that function. This time period applies to both
        the above start and stop times and hence applies twice as follows:
        * The first time period "steps" up to this event from the underlying
          function. Its mid-point is centered at `time_start`.
        * The second time period "steps" down from this event to the underlying
          function. Its mid-point is centered at `time_stop`.
    '''

    # ..................{ CONCRETE                           }..................
    def __init__(
        self, start_time: float, stop_time: float, step_width: float) -> None:
        assert types.is_numeric(start_time)
        assert types.is_numeric(stop_time)
        assert types.is_numeric(step_width)

        self.start_time = start_time
        self.stop_time = stop_time
        self.step_width = step_width
