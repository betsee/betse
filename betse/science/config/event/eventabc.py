#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all timed event classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.util.type.types import type_check, NumericTypes

# ....................{ BASE                               }....................
class SimEventABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all timed event classes.

    Instances of this class parameterize an exogenous event (e.g., cell removal,
    voltage change) to be triggered at some time step(s) of the simulation.
    '''

    # ..................{ ABSTRACT                           }..................
    @abstractmethod
    def fire(self, sim: 'betse.science.sim.Simulator', t: NumericTypes) -> None:
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
class SimEventSpikeABC(SimEventABC):
    '''
    Abstract base class of all classes describing simulation events occurring at
    only a single time step (rather than over a range of time steps).

    Attributes
    ----------------------------
    time : NumericTypes
        Time step (s) at which to trigger this action.
    _is_fired : bool
        `True` if this action's `fire()` method has already been called at a
        previous time step of the current simulation _or_ `False` otherwise.
        Defaults to `False`.
    '''

    # ..................{ CONCRETE                           }..................
    @type_check
    def __init__(self, time: NumericTypes) -> None:
        self.time = time
        self._is_fired = False

# ....................{ PULSE                              }....................
class SimEventPulseABC(SimEventABC):
    '''
    Abstract base class of all classes describing simulation events occurring
    over a range of time steps (rather than at only a single time step).

    Attributes
    ----------------------------
    start_time : NumericTypes
        Time step (s) at which to begin triggering this event.
    stop_time : NumericTypes
        Time step (s) at which to cease triggering this event.
    step_rate : NumericTypes
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

    # ..................{ CONCRETE                           }..................
    @type_check
    def __init__(
        self,
        start_time: NumericTypes,
        stop_time: NumericTypes,
        step_rate: NumericTypes,
    ) -> None:

        self.start_time = start_time
        self.stop_time = stop_time
        self.step_rate = step_rate
