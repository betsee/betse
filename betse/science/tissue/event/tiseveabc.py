#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all timed event classes.
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta  #, abstractmethod
from betse.exceptions import BetseSimEventException
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ SUPERCLASSES                      }....................
class SimEventABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all timed event classes.

    Each instance of this class parameterizes some user-defined event (e.g.,
    cell removal) to be applied at some simulation time step(s).

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
    def fire(self, phase: SimPhase, time_step: float) -> None:
        '''
        Apply this event for the passed time step of the passed simulation
        phase.

        This method notes this event to have been applied, ensuring the
        :meth:`is_fired` property to subsequently report ``True``.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        time_step : float
            Current time step of this phase being simulated.
        '''

        self._is_fired = True

# ....................{ SUBCLASSES                        }....................
class SimEventSpikeABC(SimEventABC):
    '''
    Abstract base class of all **simulation spike event** (i.e., occurring at
    only a single time step rather than over a range of time steps) subclasses.

    Attributes
    ----------
    _time_step : float
        Time step in seconds (s) at which to trigger this action.
    '''

    # ..................{ INITIALIZERS                     }..................
    @type_check
    def __init__(
        self,
        p: 'betse.science.parameters.Parameters',
        time_step: float,
    ) -> None:
        '''
        Initialize this simulation spike event for the passed simulation
        configuration.

        Parameters
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        time_step : float
            Time step in seconds (s) at which to trigger this action.

        Raises
        ----------
        BetseSimEventException
            If this time step is invalid (i.e., *not* in the range
            ``[0, p.sim_time_total)``).
        '''

        # Initialize our superclass.
        super().__init__()

        # If this time step is invalid, log a non-fatal warning.
        #
        # While an invalid time step is arguably questionable, this
        # invalidity is *NOT* fatal for most use cases and sublasses.
        # Ergo, logging a non-fatal warning is saner than raising a
        # fatal exception here.
        if not 0.0 <= time_step < p.sim_time_total:
            logs.log_warning(
                'Event time %f invalid '
                '(i.e., not in range [0.0, %f)).',
                time_step, p.sim_time_total)

        # Classify all passed parameters.
        self._time_step = time_step


class SimEventPulseABC(SimEventABC):
    '''
    Abstract base class of all **simulation pulse event** (i.e., occurring over
    a range of time steps rather than at only a single time step) subclasses.

    Attributes
    ----------
    start_time_step : float
        Time step (s) at which to begin triggering this event.
    stop_time_step : float
        Time step (s) at which to cease triggering this event.
    time_step_rate : float
        Slope of the pair of step functions guaranteeing smooth continuity
        between the background function and this event. Each step function is
        the mirror image of the other reflected across the Y axis. These are:

        * A "step" up from the background function to this event, whose
          mid-point is centered at ``start_time``.
        * A "step" down from this event back to the background function, whose
          mid-point is centered at ``stop_time``.

        If the background function is time-dependent, this slope is a **rate**
        (i.e., change over time). For the :class:`SimEventPulseVoltage`
        subclass, for example, this is the rate in voltage per seconds (V/s) at
        which:

        * The background voltage is first increased to the peak voltage.
        * The peak voltage is later decreased to the background voltage.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,
        p: 'betse.science.parameters.Parameters',
        start_time_step: float,
        stop_time_step: float,
        time_step_rate: float,
    ) -> None:
        '''
        Initialize this simulation pulse event for the passed simulation
        configuration.

        Parameters
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        start_time_step : float
            Time step (s) at which to begin triggering this event.
        stop_time_step : float
            Time step (s) at which to cease triggering this event.
        time_step_rate : float
            Slope of the pair of step functions guaranteeing smooth continuity
            between the background function and this event. See the class
            docstring for details.
        '''

        # Initialize our superclass.
        super().__init__()

        # If this start exceeds this stop time step, raise an exception.
        if start_time_step > stop_time_step:
            raise BetseSimEventException(
                'Start event time {} exceeds stop event time {}.'.format(
                    start_time_step, stop_time_step))

        # FIXME -- these might be warnings, but an exception is not feasible as these are somtimes
        # desired states.

        # # If this start or stop time step are invalid, raise an exception.
        # if not 0.0 <= start_time_step < p.sim_time_total:
        #     raise BetseSimEventException(
        #         'Start event time {} invalid (i.e., not in range '
        #         '[0.0, {})).'.format(start_time_step, p.sim_time_total))
        # if not 0.0 <= stop_time_step < p.sim_time_total:
        #     raise BetseSimEventException(
        #         'Stop event time {} invalid (i.e., not in range '
        #         '[0.0, {})).'.format(stop_time_step, p.sim_time_total))

        # Classify all passed parameters.
        self.start_time_step = start_time_step
        self.stop_time_step = stop_time_step
        self.time_step_rate = time_step_rate
