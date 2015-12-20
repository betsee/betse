#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractstaticmethod
from betse.exceptions import BetseExceptionParameters
from betse.util.io import loggers
from betse.util.type import types

# ....................{ BASE                               }....................
class Event(object, metaclass = ABCMeta):
    '''
    Abstract base class of all simulation event classes.

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
        assert types.is_numeric(start_time), \
            types.assert_not_numeric(start_time)
        assert types.is_numeric(stop_time), \
            types.assert_not_numeric(stop_time)
        assert types.is_numeric(step_width), \
            types.assert_not_numeric(step_width)

        self.start_time = start_time
        self.stop_time = stop_time
        self.step_width = step_width

# ....................{ CUT                                }....................
class EventVoltage(Event):
    '''
    Event applying a directed voltage to the environmental boundary during some
    time period of the simulation.

    A positive voltage will be applied to one boundary edge during this period;
    a negative voltage will be applied to another.

    Attributes
    ----------------------------
    peak_voltage : float
        Maximum voltage (V) to be applied.
    boundary_positive_voltage : str
        String identifying the boundary edge to apply this positive voltage to.
        Valid values include:
        * `T`, the top boundary.
        * `B`, the bottom boundary.
        * `L`, the left boundary.
        * `R`, the right boundary.
    boundary_negative_voltage : str
        String identifying the boundary edge to apply this negative voltage to.
        Valid values are as for `boundary_positive_voltage` above.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @staticmethod
    def make(params: 'Parameters') -> 'EventVoltage':
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # Object to be returned, defaulting to nothing.
        event = None

        # If this event is enabled, create an instance of this class.
        aev = params.config['apply external voltage']
        if bool(aev['event happens']):
            # If extracellular spaces are enabled, parse this event.
            if params.sim_ECM:
                event = EventVoltage(
                    start_time=float(aev['change start']),
                    stop_time=float(aev['change finish']),
                    step_width=float(aev['change rate']),
                    peak_voltage=float(aev['peak value']),
                    positive_voltage_boundary=\
                        EventVoltage._convert_boundary_str_to_char(
                            aev['boundary positive voltage']),
                    negative_voltage_boundary=\
                        EventVoltage._convert_boundary_str_to_char(
                            aev['boundary negative voltage']),
                )
            # Else, print a non-fatal warning.
            else:
                loggers.log_warning(
                    'Ignoring voltage event, '
                    'as extracellular spaces are disabled.')

        return event

        # tpd = params.config['tissue profile definition']

        # # If cut profiles are enabled, return an instance of this class.
        # if not tpd['profiles enabled']:
        #     cp = tpd['cut profile']
        #     picker = TissuePicker.make(cp['cell targets'], params)
        #     return CutProfile(cp['name'], picker)
        # # Else, return the empty void of space.
        # else:
        #     return None
        #

    # ..................{ PRIVATE ~ static                   }..................
    #FIXME: Efficiency is probably *NOT* a concern here. Ideally, BETSE
    #should use human-readable strings (e.g., "top") or perhaps even
    #enumeration constants rather than machine-readable characters
    #(e.g., "T") everywhere. Until utopia happens, this utility
    #function remains. Life to the livid givers!
    @staticmethod
    def _convert_boundary_str_to_char(boundary: str) -> str:
        '''
        Convert the passed human-readable string identifying an
        environmental boundary edge (e.g., `top`) into the
        corresponding machine-readable character (e.g., `T`).
        '''
        if boundary == 'top':
            return 'T'
        elif boundary == 'bottom':
            return 'B'
        elif boundary == 'left':
            return 'L'
        elif boundary == 'right':
            return 'R'
        else:
            raise BetseExceptionParameters(
                'Boundary edge "{}" unrecognized.'.format(boundary))

    # ..................{ PUBLIC                             }..................
    def __init__(
        self,
        start_time: float,
        stop_time: float,
        step_width: float,
        peak_voltage: float,
    ) -> None:
        assert types.is_numeric(peak_voltage), \
            types.assert_not_numeric(peak_voltage)

        super().__init__(start_time, stop_time, step_width)

        self.peak_voltage = peak_voltage
