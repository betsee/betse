#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all timed event classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod, abstractstaticmethod
from betse.util.type import types

# ....................{ BASE                               }....................
class Event(object, metaclass=ABCMeta):
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
        from the passed simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
            Current simulation configuration.

        Returns
        ----------------------------
        Event
            Concrete instance of this abstract base class.
        '''
        pass

    # ..................{ ABSTRACT                           }..................
    @abstractmethod
    def fire(self, sim: 'Simulation', t: float) -> None:
        '''
        Apply this event to the passed time step of the passed tissue
        simulation.

        Parameters
        ----------------------------
        sim : Simulation
            Current tissue simulation.
        t : float
            Time step to apply this event to.
        '''
        pass

# ....................{ PERIOD                             }....................
class Action(Event):
    '''
    Abstract base class of all classes describing simulation events occurring at
    only a single time step (rather than over a range of time steps).

    Attributes
    ----------------------------
    time : float
        Time step (s) at which to trigger this action.
    _is_fired : bool
        `True` if this action's `fire()` method has already been called at a
        previous time step of the current simulation _or_ `False` otherwise.
        Defaults to `False`.
    '''

    # ..................{ CONCRETE                           }..................
    def __init__(self, time: float) -> None:
        assert types.is_numeric(time)
        self.time = time
        self._is_fired = False

# ....................{ PULSE                              }....................
class Pulse(Event):
    '''
    Abstract base class of all classes describing simulation events occurring
    over a range of time steps (rather than at only a single time step).

    Attributes
    ----------------------------
    start_time : float
        Time step (s) at which to begin triggering this event.
    stop_time : float
        Time step (s) at which to cease triggering this event.
    step_rate : float
        Slope of the pair of step functions guaranteeing smooth continuity
        between the background function and this event. Each step function is
        the mirror image of the other reflected across the Y axis. These are:
        * A "step" up from the background function to this event, whose
          mid-point is centered at `start_time`.
        * A "step" down from this event back to the background function, whose
          mid-point is centered at `stop_time`.
        If the background function is time-dependent, this slope is a *rate*
        (i.e., change over time). For the `PulseVoltage` subclass, for
        example, this is the rate in voltage per seconds (V/s) at which:
        * The background voltage is first increased to the peak voltage.
        * The peak voltage is later decreased to the background voltage.
    '''

    # ..................{ CONCRETE                           }..................
    def __init__(
        self, start_time: float, stop_time: float, step_rate: float) -> None:
        assert types.is_numeric(start_time)
        assert types.is_numeric(stop_time)
        assert types.is_numeric(step_rate)

        self.start_time = start_time
        self.stop_time = stop_time
        self.step_rate = step_rate
